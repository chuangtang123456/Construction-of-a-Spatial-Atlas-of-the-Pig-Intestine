# clean per-sample objects -> one merged object per treatment
source("00_config.R")
suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(harmony)
})

meta  <- read.csv(METADATA_FILE, stringsAsFactors = FALSE)
group <- unique(meta[, c("Sample", "Treatment")])

merge_group <- function(g) {
  ids  <- group$Sample[group$Treatment == g]
  objs <- lapply(ids, function(i) readRDS(file.path(OUT_DIR, "clean", paste0(i, ".rds"))))

  seu <- merge(objs[[1]], objs[-1], add.cell.ids = ids)
  seu <- NormalizeData(seu, verbose = FALSE) %>%
    FindVariableFeatures(nfeatures = N_HVG, verbose = FALSE) %>%
    ScaleData(vars.to.regress = "nCount_Spatial", verbose = FALSE) %>%
    RunPCA(verbose = FALSE)
  seu <- RunHarmony(seu, group.by.vars = "sample", lambda = HARMONY_LAMBDA)
  seu <- RunUMAP(seu, dims = 1:N_PCS, reduction = "harmony")

  dir.create(file.path(OUT_DIR, "merge"), recursive = TRUE, showWarnings = FALSE)
  saveRDS(seu, file.path(OUT_DIR, "merge", paste0(g, ".rds")))
}

for (g in unique(group$Treatment)) merge_group(g)
