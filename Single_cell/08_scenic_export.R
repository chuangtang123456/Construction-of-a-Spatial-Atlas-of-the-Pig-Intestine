# export the expression matrix for pySCENIC (run before the python steps)
source("00_config.R")
suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
})

dir.create(SCENIC_DIR, recursive = TRUE, showWarnings = FALSE)

# pySCENIC runs on jejunum and colon cells together
seu <- merge(
  readRDS(file.path(OUT_DIR, "merge", "SPF_GF_FMT_Jej.rds")),
  readRDS(file.path(OUT_DIR, "merge", "SPF_GF_FMT_Col.rds"))
)
seu <- JoinLayers(seu)
seu$group_label <- paste0(seu$segment, "-", seu$cell_type)

DefaultAssay(seu) <- "RNA"
seu <- DietSeurat(seu, counts = TRUE, data = TRUE, scale.data = FALSE, assays = "RNA")

saveRDS(seu, file.path(SCENIC_DIR, "seurat_scenic.rds"))

loom <- as.loom(seu, filename = file.path(SCENIC_DIR, "input.loom"), verbose = FALSE)
loom$close_all()
