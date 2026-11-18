# merge samples into the four analysis groups and integrate each one
source("00_config.R")
suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(harmony)
})

merge_group <- function(name) {
  cfg <- MERGE_GROUPS[[name]]

  ids  <- SAMPLE_MAP$sample[SAMPLE_MAP$segment == cfg$segment &
                              SAMPLE_MAP$treatment %in% cfg$treatment]
  objs <- lapply(ids, function(i) readRDS(file.path(OUT_DIR, "clean", paste0(i, ".rds"))))

  seu <- merge(objs[[1]], objs[-1], add.cell.ids = ids)
  seu <- NormalizeData(seu, verbose = FALSE) %>%
    FindVariableFeatures(nfeatures = N_HVG, verbose = FALSE) %>%
    ScaleData(vars.to.regress = "nCount_RNA", verbose = FALSE) %>%
    RunPCA(verbose = FALSE)
  seu <- RunHarmony(seu, group.by.vars = "sample", lambda = cfg$lambda)
  seu <- RunUMAP(seu, dims = 1:cfg$dims, reduction = "harmony")
  seu <- FindNeighbors(seu, reduction = "harmony", dims = 1:cfg$dims) %>%
    FindClusters(resolution = cfg$resolution)

  dir.create(file.path(OUT_DIR, "merge"), recursive = TRUE, showWarnings = FALSE)
  saveRDS(seu, file.path(OUT_DIR, "merge", paste0(name, ".rds")))
}

for (name in names(MERGE_GROUPS)) merge_group(name)
