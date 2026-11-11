# RCTD deconvolution: spot cell type proportions from a single-cell reference
source("00_config.R")
suppressPackageStartupMessages({
  library(Seurat)
  library(spacexr)
})

refs <- readRDS(SC_REF_FILE)
meta <- read.csv(METADATA_FILE, stringsAsFactors = FALSE)

make_reference <- function(ref) {
  keep <- names(which(table(ref$Cell.type) >= RCTD_MIN_CELLS))
  ref  <- subset(ref, cells = colnames(ref)[ref$Cell.type %in% keep])

  counts <- round(GetAssayData(ref, layer = "counts"))
  ct <- factor(ref$Cell.type)
  names(ct) <- colnames(counts)

  Reference(counts, ct, colSums(counts))
}

ref_jej <- make_reference(refs$jej)
ref_col <- make_reference(refs$col)

rctd_one <- function(id) {
  seu <- readRDS(file.path(OUT_DIR, "clean", paste0(id, ".rds")))
  ref <- if (unique(seu$segment) %in% c("Duo", "Jej", "Ile")) ref_jej else ref_col

  seu <- subset(seu, cells = colnames(seu)[seu$layer %in% RCTD_LAYERS])

  counts <- GetAssayData(seu, layer = "counts")
  xy     <- GetTissueCoordinates(seu, cols = c("imagerow", "imagecol"))[colnames(seu), ]
  coords <- data.frame(x = xy[, 1], y = xy[, 2], row.names = colnames(seu))

  puck <- SpatialRNA(coords, counts, colSums(counts))
  rctd <- create.RCTD(puck, ref, max_cores = RCTD_MAX_CORES)
  rctd <- run.RCTD(rctd, doublet_mode = "full")

  w <- t(rctd@results$weights)
  colnames(w) <- sub("^weights\\.", "", colnames(w))

  dir.create(file.path(OUT_DIR, "rctd"), recursive = TRUE, showWarnings = FALSE)
  write.csv(w, file.path(OUT_DIR, "rctd", paste0(id, ".csv")))
}

for (id in unique(meta$Sample)) rctd_one(id)
