# LTSR: cell type composition across treatments, per segment
source("00_config.R")
suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(ggplot2)
  library(lme4)
  library(numDeriv)
})
source("R/ltsr_functions.R")

out_dir <- file.path(OUT_DIR, "ltsr")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

ltsr_one <- function(segment) {
  seu <- readRDS(file.path(OUT_DIR, "merge", paste0("SPF_GF_FMT_", segment, ".rds")))

  obs <- data.frame(
    Sample_ID = seu$sample,
    Treatment = seu$treatment,
    Cell_type = seu$cell_type
  )

  res <- CellTypeCompositionAnalysis(
    obs,
    colSample = "Sample_ID", colCelltype = "Cell_type", colVarCats = "Treatment"
  )

  write.csv(res$ranef, file.path(out_dir, paste0(segment, "_ranef.csv")), row.names = FALSE)
  write.csv(res$sdse,  file.path(out_dir, paste0(segment, "_sdse.csv")),  row.names = FALSE)

  p_ranef <- plot_ranef(
    res$ranef,
    vars = list(Treatment = LTSR_GROUPS),
    celltypes = sort(unique(obs$Cell_type)),
    references = list(Treatment = LTSR_REFERENCE),
    celltype_order = "hclust",
    maxFC = LTSR_MAXFC, LTSR2p = FALSE
  )
  p_sdse <- plot_sdse(res$sdse, "Sample_ID")

  ggsave(file.path(out_dir, paste0(segment, "_ranef.pdf")), p_ranef, width = 6, height = 7)
  ggsave(file.path(out_dir, paste0(segment, "_sdse.pdf")),  p_sdse,  width = 4, height = 3)
}

for (segment in LTSR_SEGMENTS) ltsr_one(segment)
