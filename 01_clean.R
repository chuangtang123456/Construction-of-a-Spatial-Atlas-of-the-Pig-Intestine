# raw spaceranger output + metadata table -> clean per-sample object
source("00_config.R")
suppressPackageStartupMessages({
  library(Seurat)
})

meta <- read.csv(METADATA_FILE, stringsAsFactors = FALSE)

clean_one <- function(id) {
  seu <- Load10X_Spatial(file.path(RAW_DIR, id))
  seu[["percent.mt"]] <- PercentageFeatureSet(seu, features = MT_GENES)
  seu <- subset(seu, subset = nCount_Spatial > UMI_MIN & percent.mt < MT_MAX)

  counts <- GetAssayData(seu, layer = "counts")
  seu <- subset(seu, features = rownames(counts)[rowSums(counts > 1) >= MIN_SPOTS_PER_GENE])

  m <- meta[meta$Sample == id, ]
  seu <- subset(seu, cells = intersect(colnames(seu), m$Cell_id))
  m <- m[match(colnames(seu), m$Cell_id), ]

  seu$layer     <- m$Layer
  seu$layer2    <- m$Layer2
  seu$zone      <- m$Zone
  seu$unroll_x  <- m$Unroll_x
  seu$unroll_y  <- m$Unroll_y
  seu$sample    <- id
  seu$treatment <- m$Treatment
  seu$segment   <- m$Segment

  dir.create(file.path(OUT_DIR, "clean"), recursive = TRUE, showWarnings = FALSE)
  saveRDS(seu, file.path(OUT_DIR, "clean", paste0(id, ".rds")))
}

for (id in unique(meta$Sample)) clean_one(id)
