# per sample filtering: ambient RNA (SoupX), doublets (scDblFinder), cell QC,
# then attach the given annotation
source("00_config.R")
suppressPackageStartupMessages({
  library(Seurat)
  library(SoupX)
  library(scDblFinder)
  library(SingleCellExperiment)
  library(dplyr)
})

meta <- bind_rows(lapply(METADATA_FILES, read.csv, stringsAsFactors = FALSE))

filter_one <- function(folder) {
  sample <- SAMPLE_MAP$sample[SAMPLE_MAP$folder == folder]

  # ambient RNA ------------------------------------------------------------
  raw <- Read10X(file.path(RAW_DIR, folder, "raw"))
  fil <- Read10X(file.path(RAW_DIR, folder, "filtered"))

  seu <- CreateSeuratObject(fil, project = folder) %>%
    NormalizeData(verbose = FALSE) %>%
    FindVariableFeatures(verbose = FALSE) %>%
    ScaleData(verbose = FALSE) %>%
    RunPCA(verbose = FALSE) %>%
    FindNeighbors(dims = 1:N_PCS, verbose = FALSE) %>%
    FindClusters(resolution = SOUP_RESOLUTION, verbose = FALSE)

  ch <- SoupChannel(raw, fil, calcSoupProfile = FALSE) %>%
    estimateSoup() %>%
    setClusters(clusters = seu$seurat_clusters) %>%
    autoEstCont(doPlot = FALSE)

  seu <- CreateSeuratObject(counts = adjustCounts(ch)) %>%
    NormalizeData(verbose = FALSE) %>%
    FindVariableFeatures(verbose = FALSE) %>%
    ScaleData(verbose = FALSE) %>%
    RunPCA(verbose = FALSE) %>%
    FindNeighbors(dims = 1:N_PCS, verbose = FALSE) %>%
    FindClusters(resolution = SOUP_RESOLUTION, verbose = FALSE)

  # doublets ---------------------------------------------------------------
  sce <- scDblFinder(as.SingleCellExperiment(seu), clusters = "seurat_clusters")
  seu$doublet <- sce$scDblFinder.class
  seu <- subset(seu, cells = colnames(seu)[seu$doublet == "singlet"])

  # cell QC ----------------------------------------------------------------
  seu[["percent.mt"]] <- PercentageFeatureSet(seu, features = MT_GENES)
  seu <- subset(seu, subset = nCount_RNA > UMI_MIN & percent.mt < MT_MAX)

  # annotation -------------------------------------------------------------
  m <- meta[meta$Sample == sample, ]
  seu <- subset(seu, cells = intersect(colnames(seu), m$Cell_id))
  m <- m[match(colnames(seu), m$Cell_id), ]

  seu$cell_type <- m$Cell.type
  seu$sample    <- sample
  seu$treatment <- m$Treatment
  seu$segment   <- m$Segment
  seu$group     <- m$Group

  seu <- DietSeurat(seu, layers = "counts")

  dir.create(file.path(OUT_DIR, "clean"), recursive = TRUE, showWarnings = FALSE)
  saveRDS(seu, file.path(OUT_DIR, "clean", paste0(sample, ".rds")))
}

for (folder in SAMPLE_MAP$folder) filter_one(folder)
