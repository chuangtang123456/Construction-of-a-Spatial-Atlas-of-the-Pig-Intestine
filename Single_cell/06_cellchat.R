# CellChat: cell-cell communication per treatment and segment, then comparisons
source("00_config.R")
suppressPackageStartupMessages({
  library(Seurat)
  library(CellChat)
  library(dplyr)
  library(ggplot2)
  library(patchwork)
})

out_dir <- file.path(OUT_DIR, "cellchat")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# pig symbols -> human symbols, one gene per symbol
homo <- read.csv(ORTHOLOG_FILE, check.names = FALSE, stringsAsFactors = FALSE)
homo <- homo[homo[["Human homology type"]] == "ortholog_one2one" & !is.na(homo[["Human gene name"]]), ]
human_of <- tapply(homo[["Human gene name"]], homo[["Gene name"]], function(x) x[1])

build_one <- function(segment, treatment) {
  seu <- readRDS(file.path(OUT_DIR, "merge", paste0("SPF_GF_FMT_", segment, ".rds")))
  seu <- subset(seu, cells = colnames(seu)[seu$treatment == treatment])

  mat <- GetAssayData(seu, layer = "counts")
  mat <- mat[!is.na(human_of[rownames(mat)]), ]
  mat <- mat[order(rowSums(mat), decreasing = TRUE), ]
  rownames(mat) <- human_of[rownames(mat)]
  mat <- mat[!duplicated(rownames(mat)), ]

  seu <- CreateSeuratObject(counts = mat, meta.data = seu@meta.data, assay = "RNA")
  seu <- NormalizeData(seu, verbose = FALSE) %>%
    FindVariableFeatures(verbose = FALSE) %>%
    ScaleData(verbose = FALSE)

  meta <- data.frame(row.names = colnames(seu), label = seu$cell_type, samples = seu$sample)

  cellchat <- createCellChat(object = GetAssayData(seu, layer = "data"), meta = meta, group.by = "label")
  cellchat@DB <- subsetDB(CellChatDB.human, search = CELLCHAT_DB, key = "annotation")

  cellchat <- subsetData(cellchat)
  cellchat <- identifyOverExpressedGenes(cellchat)
  cellchat <- identifyOverExpressedInteractions(cellchat, variable.both = FALSE)

  cellchat <- computeCommunProb(
    cellchat,
    type = "truncatedMean", trim = 0.1,
    distance.use = TRUE, interaction.range = 250, scale.distance = 0.01,
    contact.dependent = TRUE, contact.range = 100
  )
  cellchat <- filterCommunication(cellchat, min.cells = CELLCHAT_MIN_CELLS)
  cellchat <- computeCommunProbPathway(cellchat)
  cellchat <- aggregateNet(cellchat)

  saveRDS(cellchat, file.path(out_dir, paste0(treatment, "_", segment, ".rds")))
}

compare_two <- function(segment, comparison) {
  groups <- CELLCHAT_COMPARISONS[[comparison]]
  objects <- lapply(groups, function(g) readRDS(file.path(out_dir, paste0(g, "_", segment, ".rds"))))
  names(objects) <- tolower(groups)

  merged <- mergeCellChat(objects, add.names = names(objects))

  p_count <- compareInteractions(merged, show.legend = FALSE, group = names(objects)) +
    compareInteractions(merged, show.legend = FALSE, group = names(objects), measure = "weight")
  p_heat <- netVisual_heatmap(merged) + netVisual_heatmap(merged, measure = "weight")
  p_rank <- rankNet(merged, mode = "comparison", stacked = TRUE, do.stat = TRUE)

  tag <- paste0(segment, "_", comparison)
  ggsave(file.path(out_dir, paste0(tag, "_counts.pdf")), p_count, width = 8, height = 4)
  ggsave(file.path(out_dir, paste0(tag, "_heatmap.pdf")), p_heat, width = 12, height = 5)
  ggsave(file.path(out_dir, paste0(tag, "_ranknet.pdf")), p_rank, width = 5, height = 5)
}

future::plan("multisession", workers = CELLCHAT_WORKERS)

for (segment in CELLCHAT_SEGMENTS) {
  for (treatment in CELLCHAT_TREATMENTS) build_one(segment, treatment)
  for (comparison in names(CELLCHAT_COMPARISONS)) compare_two(segment, comparison)
}
