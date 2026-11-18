# Single_cell

单细胞：原始数据 → 过滤 → 合并成四个分析对象 → 组成 / 邻域丰度 / 差异基因 / 通讯 / 注释可重复性 / 调控网络。

| 脚本 | 输入 | 输出 |
|---|---|---|
| `00_config.R` | — | 路径、样本映射、QC、整合、LTSR、milo、DEG、CellChat、MetaNeighbor 与 SCENIC 参数 |
| `01_filter.R` | `RAW_DIR/<folder>/{raw,filtered}` + metadata | `output/clean/<sample>.rds` |
| `02_merge.R` | `output/clean/*.rds` | `output/merge/<group>.rds` |
| `03_ltsr.R` | `output/merge/SPF_GF_FMT_*.rds` | `output/ltsr/<segment>_{ranef,sdse}.csv`、`.pdf` |
| `04_milo.R` | `output/merge/SPF_GF_FMT_*.rds` | `output/milo/<segment>_<cmp>_da.csv`、`_graph.pdf`、`_sina.pdf` |
| `05_deg.R` | `output/merge/SPF_GF_FMT_*.rds` | `output/deg/<segment>_{all,clean}.csv`、`_counts.pdf` |
| `06_cellchat.R` | `output/merge/SPF_GF_FMT_*.rds` | `output/cellchat/<treatment>_<segment>.rds`、比较图 |
| `07_metaneighbor.R` | `output/merge/SPF_GF_FMT_*.rds` | `output/metaneighbor/<tag>_auroc.csv`、`.pdf` |
| `08_scenic_export.R` | `output/merge/SPF_GF_FMT_*.rds` | `output/scenic/input.loom`、`seurat_scenic.rds` |
| `09_scenic_run.sh` | `input.loom` + TF 列表 | `adj.tsv`、`reg.csv`、`scenic.loom` |
| `10_scenic_import.R` | pySCENIC 输出的 loom | `output/scenic/rss.csv`、`rss_top.pdf`、`rss_heatmap.pdf` |
| `R/ltsr_functions.R` | — | LTSR 方法函数（`CellTypeCompositionAnalysis`、`plot_ranef`、`plot_sdse`） |

## 样本名映射

磁盘目录名和 metadata 里的 `Sample` 不一致，`00_config.R` 的 `SAMPLE_MAP` 负责对应：

| 目录 | Sample |
|---|---|
| `SPF-A-K-1` … `SPF-C-K-1` | `SPF_Jej_A_1` … `SPF_Jej_C_1` |
| `SPF-A-J-1` … `SPF-C-J-1` | `SPF_Col_A_1` … `SPF_Col_C_1` |
| `A_KST_1` … `C_KST_1` | `GF_Jej_A_1` … `GF_Jej_C_1` |
| `A_JST_1` … `C_JST_1` | `GF_Col_A_1` … `GF_Col_C_1` |
| `FMT-KSC-A1/B1/C1` | `FMT_Jej_A_1/B_1/C_1` |
| `FMT-JSC-A1`, `FMT-JSC-C1`, `FMT-JST-B1` | `FMT_Col_A_1`, `FMT_Col_C_1`, `FMT_Col_B_1` |

`treatment` 和 `segment` 两列从 `sample` 名解析出来（`SPF_Jej_A_1` → SPF / Jej）。

## 过滤（`01_filter.R`）

一个函数跑完四件事：

1. **SoupX**：`raw/` 与 `filtered/` 建 `SoupChannel`，粗聚类（`SOUP_RESOLUTION = 0.5`）估环境 RNA，`adjustCounts()` 去污染；
2. **scDblFinder**：去污染后重聚类并判双细胞，只留 `singlet`；
3. **QC**：`nCount_RNA > 2000`、`percent.mt < 35`；
4. **注释**：按 `Cell_id` 从 `Metadata_SC_*.csv` 取 `Cell.type`，以及 `Treatment` / `Segment` / `Group`。

最后 `DietSeurat(layers = "counts")` 只留 counts，避免 18 份 scale.data 占内存。

## 合并（`02_merge.R`）

四组分别合并、分别整合，参数在 `00_config.R` 的 `MERGE_GROUPS` 里逐组写死：

| 组 | 样本 | harmony `lambda` | UMAP `dims` | `resolution` |
|---|---|---:|---:|---:|
| `SPF_Jej` | 3 个 SPF 空肠 | 2 | 20 | 0.1 |
| `SPF_Col` | 3 个 SPF 结肠 | 2 | 20 | 0.1 |
| `SPF_GF_FMT_Jej` | 9 个空肠 | 2 | 30 | 0.1 |
| `SPF_GF_FMT_Col` | 9 个结肠 | 2 | 30 | 0.1 |

每组流程一致：`NormalizeData` → HVG 2000 → `ScaleData(nCount_RNA)` → PCA → harmony → UMAP → `FindNeighbors`/`FindClusters`。输出 `output/merge/<group>.rds`。

刻意剔除特定细胞簇的步骤（如原来结肠里去掉 `ENPP7+` 污染群）不属于标准流程，这版不包含。

`dims` 的差别来自原脚本：三处理组一起整合时用 30，SPF 单独整合时用 20。`FindNeighbors` 现在显式用 `dims = 1:dims`（原脚本没传，走的是默认 1:10）。

## LTSR 组成分析（`03_ltsr.R`）

对 `SPF_GF_FMT_Jej`、`SPF_GF_FMT_Col` 两个对象各跑一次，检验细胞类型比例在 SPF / GF / FMT 之间的差异。

1. 从对象里取 `sample`、`treatment`、`cell_type` 三列，构成每个细胞的观测表；
2. `CellTypeCompositionAnalysis()` 拟 GLMM：`Treatment` 为固定效应，样本、细胞类型及其嵌套为随机效应，泊松分布；
3. `plot_ranef()` 画随机效应 dotplot（以 `GF` 为参照，`maxFC = 3`），`plot_sdse()` 画方差解释的 forest plot。

| 输出 | 内容 |
|---|---|
| `output/ltsr/<segment>_ranef.csv` | 随机效应后验均值与标准差 |
| `output/ltsr/<segment>_sdse.csv` | 方差参数及其标准误 |
| `output/ltsr/<segment>_ranef.pdf` | 细胞类型 × 处理组的效应 dotplot |
| `output/ltsr/<segment>_sdse.pdf` | 各变量解释方差的 forest plot |

`R/ltsr_functions.R` 是该方法自带的函数库（498 行，来自原脚本），只保留函数部分，去掉了写 `*.RData` 的副作用；`.make_formula()` 里加了一个保护，允许只有一个分类变量（原版至少要两个，否则 `combn()` 报错）。

## miloR 邻域丰度（`04_milo.R`）

在每个肠段内做两次两两比较：`SPF vs GF`、`FMT vs GF`（GF 为参照，`logFC > 0` 表示在测试组富集）。

1. 从 `SPF_GF_FMT_<segment>` 中取出两组样本，`treatment` 设为因子（参照在前）；
2. `Milo()` → `buildGraph(k = 30, d = 30, reduced.dim = "harmony")` → `makeNhoods(prop = 0.1, refined = TRUE)` → `countCells(sample)` → `calcNhoodDistance()`；
3. `design.df` 用每个样本的 `sample × treatment`，`testNhoods(design = ~ treatment)`；
4. `annotateNhoods(coldata_col = "cell_type")` 给邻域打细胞类型标签，按 `SpatialFDR < 0.2` 和 `logFC` 方向分组。

| 输出 | 内容 |
|---|---|
| `output/milo/<segment>_<cmp>_da.csv` | 每个邻域的 logFC、SpatialFDR、所属细胞类型、分组 |
| `output/milo/<segment>_<cmp>_graph.pdf` | 邻域图，按 logFC 着色 |
| `output/milo/<segment>_<cmp>_sina.pdf` | 按细胞类型的 logFC sina 图 |

参数在 `00_config.R`：`MILO_K`、`MILO_D`、`MILO_PROP`、`MILO_REDUCTION`、`MILO_FDR`、`MILO_COMPARISONS`。

## 处理间差异基因（`05_deg.R`）

在每个肠段、每种细胞类型内做两次比较：`SPF vs GF`、`FMT vs GF`（GF 为参照，`avg_log2FC > 0` 表示在测试组更高）。

1. 细胞类型必须三个处理组都有 ≥ `DEG_MIN_CELLS = 10` 个细胞才进入检验；
2. 每个「细胞类型 × 比较」单独 `FindMarkers(ident.1 = 测试组, ident.2 = GF)`；
3. 显著阈值：`p_val_adj < 0.05`、`pct.1 > 0.3 | pct.2 > 0.3`、`|avg_log2FC| > log2(1.5)`；
4. `shared` 列标出两次比较里方向一致的基因（SPF 和 FMT 相对 GF 同向变化）。

| 输出 | 内容 |
|---|---|
| `output/deg/<segment>_all.csv` | 全部检验结果（未过滤） |
| `output/deg/<segment>_clean.csv` | 过滤后的 DEG，含 `shared` 标记 |
| `output/deg/<segment>_counts.pdf` | 每个细胞类型、每个比较的 DEG 数量柱状图 |

## 细胞通讯（`06_cellchat.R`）

对每个「肠段 × 处理」建一个 CellChat 对象，再做两两比较。

1. 从 `SPF_GF_FMT_<segment>` 取出单个处理的细胞，基因名换成人类符号（CellChatDB 是人源库）；
2. `createCellChat(group.by = "label")` → 用 `CellChatDB.human` 的 `Secreted Signaling` 子库 → `subsetData()` → `identifyOverExpressedGenes()` → `identifyOverExpressedInteractions()`；
3. `computeCommunProb(type = "truncatedMean", trim = 0.1, ...)` → `filterCommunication(min.cells = 3)` → `computeCommunProbPathway()` → `aggregateNet()`；
4. 比较：`mergeCellChat()` 后跑 `compareInteractions()`、`netVisual_heatmap()`、`rankNet()`。

| 输出 | 内容 |
|---|---|
| `output/cellchat/<treatment>_<segment>.rds` | 单处理 CellChat 对象（6 个） |
| `output/cellchat/<segment>_<cmp>_counts.pdf` | 交互数量 / 强度对比 |
| `output/cellchat/<segment>_<cmp>_heatmap.pdf` | 细胞类型间通讯热图（数量 + 强度） |
| `output/cellchat/<segment>_<cmp>_ranknet.pdf` | 通路层面的差异通路排序 |

参数：`CELLCHAT_DB`、`CELLCHAT_MIN_CELLS`、`CELLCHAT_WORKERS`、`CELLCHAT_COMPARISONS`。`computeCommunProb` 里的 `distance.use` / `interaction.range` 沿用原脚本；scRNA 对象没有空间坐标时 CellChat 会自动忽略距离相关参数。

## 注释可重复性（`07_metaneighbor.R`）

用无监督 `MetaNeighborUS` 检验同一套细胞类型注释在不同条件之间能否互相复现：给定 `study_id`（分组）和 `cell_type`，得到 cell type × cell type 的 AUROC 矩阵。

| 分组 | 对象 | `study_id` | 输出前缀 |
|---|---|---|---|
| 处理之间 | `SPF_GF_FMT_Jej`、`SPF_GF_FMT_Col` 各一次 | `treatment` | `<segment>_treatment` |
| 肠段之间 | 两个对象合并 | `segment` | `segment` |

流程：`variableGenes(exp_labels = study_id)` 挑基因 → `MetaNeighborUS(fast_version = TRUE)` → AUROC 矩阵 → `plotHeatmap()`。

| 输出 | 内容 |
|---|---|
| `output/metaneighbor/<tag>_auroc.csv` | AUROC 相似度矩阵（行/列名 `分组|细胞类型`） |
| `output/metaneighbor/<tag>_auroc.pdf` | 热图，看对角块是否成块 |

参数：`MN_SEGMENTS`。原脚本里还有一版用 marker 基因表代替 `variableGenes()`，没搬；需要的话把 `var_genes` 换成读表即可。

## 调控网络（`08_scenic_export.R` → `09_scenic_run.sh` → `10_scenic_import.R`）

SCENIC 分三段，中间那步是 shell，必须在装了 pyscenic 的 Python 环境里跑。

**导出**：把两个肠段的对象合并，标记 `group_label = Segment-CellType`，`DietSeurat(counts/data, scale.data = FALSE)` 后 `as.loom()` 写出 `input.loom`，同时存 `seurat_scenic.rds` 给最后一步对齐元数据。

**pySCENIC**：`09_scenic_run.sh`，三步对应原来的 `02_scenic.sh`：

```bash
pyscenic grn --num_workers 10 --output adj.tsv --method grnboost2 \
  input.loom hs_hgnc_curated_tfs2pig.txt

pyscenic ctx adj.tsv <motif rankings.feather> \
  --annotations_fname motifs-v9-nr.hgnc-m0.001-o0.0.tbl \
  --expression_mtx_fname input.loom --mode dask_multiprocessing \
  --output reg.csv --num_workers 10 --mask_dropouts

pyscenic aucell input.loom reg.csv --output scenic.loom --num_workers 10
```

脚本顶部集中了 `PROJECT_DIR`、`WORKERS = 10`、TF 列表路径，以及两个 Aertslab 资源文件名（`hg38__refseq-r80__10kb_up_and_down_tss.mc9nr.genes_vs_motifs.rankings.feather`、`motifs-v9-nr.hgnc-m0.001-o0.0.tbl`），下载地址写在注释里。

**导入**：读 `scenic.loom` 的 `RegulonsAUC` → `getAUC()` 得到「regulon × 细胞」→ `SCENIC::calcRSS(AUC, cellAnnotation = group_label)`。

| 输出 | 内容 |
|---|---|
| `output/scenic/rss.csv` | regulon × 分组 的 RSS 矩阵 |
| `output/scenic/rss_top.pdf` | 每个分组 RSS 最高的 10 个 regulon |
| `output/scenic/rss_heatmap.pdf` | regulon 之间的 RSS 相关性热图 |

参数：`SCENIC_DIR`、`SCENIC_LOOM`（pySCENIC 输出的路径）、`SCENIC_TF_LIST`、`SCENIC_WORKERS`、`SCENIC_TOP = 10`。

## 运行

```r
setwd("E:/my_project/70_intestine_code/Single_cell")
source("01_filter.R")
source("02_merge.R")
source("03_ltsr.R")
source("04_milo.R")
source("05_deg.R")
source("06_cellchat.R")
source("07_metaneighbor.R")
source("08_scenic_export.R")   # 然后 bash 09_scenic_run.sh
source("10_scenic_import.R")
```

`01_filter.R` 对 18 个样本逐个跑 SoupX，耗时最长；建议先跑 `SAMPLE_MAP$folder[1]` 确认环境正常，再整体跑。
