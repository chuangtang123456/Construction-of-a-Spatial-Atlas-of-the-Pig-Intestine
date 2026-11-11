# paths
PROJECT_DIR <- "E:/my_project/70_intestine_code"
RAW_DIR     <- "E:/my_project/57_TC_doc/01_spatial/00_rawdata"
OUT_DIR     <- file.path(PROJECT_DIR, "spatial", "output")

# metadata: one row per spot, columns Sample / Cell_id / Layer / Layer2 / Zone /
# Unroll_x / Unroll_y / Treatment / Segment
METADATA_FILE <- file.path(PROJECT_DIR, "01_rawdata", "Metadata_ST.csv")

# qc
UMI_MIN            <- 2000
MT_MAX             <- 15
MIN_SPOTS_PER_GENE <- 10
MT_GENES <- c("ND6", "ND5", "ND4L", "ND4", "ND3", "ND2", "ND1",
              "CYTB", "COX3", "COX2", "COX1", "ATP8", "ATP6")

# integration
N_HVG          <- 2000
N_PCS          <- 20
HARMONY_LAMBDA <- 1

# differential expression
MIN_PCT      <- 0.3
MIN_DIFF_PCT <- 0.1
LOGFC_MIN    <- log2(1.5)
ONLY_POS     <- TRUE

# RCTD deconvolution
# annotated single-cell reference, list with $jej (small intestine) and
# $col (large intestine), cell type stored in Cell.type
SC_REF_FILE    <- "E:/my_project/66_TC/rawdata/sc/for.plot.SPF.GF.each.rds"
RCTD_LAYERS    <- c("Mucosa", "Submucosa")
RCTD_MIN_CELLS <- 25
RCTD_MAX_CORES <- 10

# spatial aggregation of cell types
CELL_DIST_QUANTILE <- 0.5
CELL_DIST_K        <- 6

# zone programs
# the two SPF jejunum samples drive the cross-zone gene clustering
ZONE_SAMPLES <- c("SPF02_KST_2", "SPF04_KST_1")

NMF_RANK <- 9
NMF_HVG  <- 5000

ZONE_MIN_PCT <- 0.3
ZONE_PADJ    <- 0.05
ZONE_LOG2FC  <- log2(1.5)
SPARK_PADJ   <- 0.05
SPARK_CORES  <- 1

HCLUST_HEIGHT <- 1.2

# pathway activity
# pig -> human one-to-one orthologs; genes are renamed to human symbols because
# PROGENy ships a human model
ORTHOLOG_FILE <- file.path(PROJECT_DIR, "01_rawdata", "pig2huam2mouse_ensembl106.csv")
PROGENY_TOP   <- 300

# cross-species zonation
# human / mouse villus zonation means (Buercalff 2022, Moor 2018)
CROSS_SPECIES_REF <- file.path(PROJECT_DIR, "01_rawdata", "41586_2024_7793_MOESM8_ESM.csv")
CROSS_SAMPLES     <- ZONE_SAMPLES
CROSS_MIN_EXPR    <- 5e-5
CROSS_CENTER_CUT  <- 0.25
CROSS_PEARSON_CUT <- -0.5
