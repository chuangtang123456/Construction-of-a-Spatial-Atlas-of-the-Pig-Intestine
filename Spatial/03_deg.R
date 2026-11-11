# DEGs between segments (within each layer) and between layers
source("00_config.R")
suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
})

deg_dir <- file.path(OUT_DIR, "deg")
dir.create(deg_dir, recursive = TRUE, showWarnings = FALSE)

groups <- unique(read.csv(METADATA_FILE, stringsAsFactors = FALSE)$Treatment)

find_deg <- function(seu, group_by) {
  Idents(seu) <- seu@meta.data[[group_by]]
  FindAllMarkers(
    seu,
    only.pos = ONLY_POS, logfc.threshold = LOGFC_MIN,
    min.pct = MIN_PCT, min.diff.pct = MIN_DIFF_PCT, verbose = FALSE
  )
}

for (g in groups) {
  seu <- JoinLayers(readRDS(file.path(OUT_DIR, "merge", paste0(g, ".rds"))))

  # between layers
  write.csv(find_deg(seu, "layer2"),
            file.path(deg_dir, sprintf("layer_%s.csv", g)), row.names = FALSE)

  # between segments, inside each layer
  for (lyr in na.omit(unique(seu$layer2))) {
    cells <- colnames(seu)[!is.na(seu$layer2) & seu$layer2 == lyr]
    if (length(unique(seu$segment[cells])) < 2) next

    sub <- subset(seu, cells = cells)
    write.csv(find_deg(sub, "segment"),
              file.path(deg_dir, sprintf("segment_%s_%s.csv", g, lyr)), row.names = FALSE)
  }
}
