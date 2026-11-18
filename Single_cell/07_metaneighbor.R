# MetaNeighbor: are the cell type labels reproducible across conditions?
# run twice - across treatments (within a segment) and across segments
source("00_config.R")
suppressPackageStartupMessages({
  library(Seurat)
  library(MetaNeighbor)
  library(SingleCellExperiment)
})

out_dir <- file.path(OUT_DIR, "metaneighbor")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

as_sce <- function(seu) {
  seu <- JoinLayers(seu)
  SingleCellExperiment(
    assays = list(counts = GetAssayData(seu, layer = "counts")),
    colData = seu@meta.data
  )
}

run_mn <- function(sce, study, tag) {
  var_genes <- variableGenes(dat = sce, exp_labels = sce[[study]])

  auroc <- MetaNeighborUS(
    var_genes = var_genes, dat = sce,
    study_id = sce[[study]], cell_type = sce$cell_type,
    fast_version = TRUE
  )

  write.csv(auroc, file.path(out_dir, paste0(tag, "_auroc.csv")))

  pdf(file.path(out_dir, paste0(tag, "_auroc.pdf")), width = 10, height = 10)
  print(plotHeatmap(auroc))
  dev.off()
}

# across treatments, one run per segment
for (segment in MN_SEGMENTS) {
  seu <- readRDS(file.path(OUT_DIR, "merge", paste0("SPF_GF_FMT_", segment, ".rds")))
  run_mn(as_sce(seu), "treatment", paste0(segment, "_treatment"))
}

# across segments, both segments in one object
seu <- merge(
  readRDS(file.path(OUT_DIR, "merge", "SPF_GF_FMT_Jej.rds")),
  readRDS(file.path(OUT_DIR, "merge", "SPF_GF_FMT_Col.rds"))
)
run_mn(as_sce(seu), "segment", "segment")
