# paths
PROJECT_DIR <- "E:/my_project/70_intestine_code"
RAW_DIR     <- "E:/my_project/57_TC_doc/02_scell/00rawdata"
OUT_DIR     <- file.path(PROJECT_DIR, "Single_cell", "output")

# metadata: one row per cell, split by segment
METADATA_FILES <- file.path(PROJECT_DIR, "Sample_metdata",
                            c("Metadata_SC_Jej.csv", "Metadata_SC_Col.csv"))

# folder on disk -> sample name used in the metadata
SAMPLE_MAP <- data.frame(
  folder = c("SPF-A-K-1", "SPF-B-K-1", "SPF-C-K-1",
             "SPF-A-J-1", "SPF-B-J-1", "SPF-C-J-1",
             "A_KST_1",   "B_KST_1",   "C_KST_1",
             "A_JST_1",   "B_JST_1",   "C_JST_1",
             "FMT-KSC-A1", "FMT-KSC-B1", "FMT-KSC-C1",
             "FMT-JSC-A1", "FMT-JSC-C1", "FMT-JST-B1"),
  sample = c("SPF_Jej_A_1", "SPF_Jej_B_1", "SPF_Jej_C_1",
             "SPF_Col_A_1", "SPF_Col_B_1", "SPF_Col_C_1",
             "GF_Jej_A_1",  "GF_Jej_B_1",  "GF_Jej_C_1",
             "GF_Col_A_1",  "GF_Col_B_1",  "GF_Col_C_1",
             "FMT_Jej_A_1", "FMT_Jej_B_1", "FMT_Jej_C_1",
             "FMT_Col_A_1", "FMT_Col_C_1", "FMT_Col_B_1"),
  stringsAsFactors = FALSE
)
SAMPLE_MAP$treatment <- sub("_.*$", "", SAMPLE_MAP$sample)
SAMPLE_MAP$segment   <- sub("_.*$", "", sub("^[^_]+_", "", SAMPLE_MAP$sample))

# qc
UMI_MIN  <- 2000
MT_MAX   <- 35
MT_GENES <- c("ND6", "ND5", "ND4L", "ND4", "ND3", "ND2", "ND1",
              "CYTB", "COX3", "COX2", "COX1", "ATP8", "ATP6")

# ambient RNA / doublets
SOUP_RESOLUTION <- 0.5
N_PCS           <- 30

# integration: one object per group
# treatment, segment, harmony lambda, UMAP dimensions, clustering resolution
N_HVG <- 2000

MERGE_GROUPS <- list(
  SPF_Jej = list(
    treatment = "SPF", segment = "Jej",
    lambda = 2, dims = 20, resolution = 0.1
  ),
  SPF_Col = list(
    treatment = "SPF", segment = "Col",
    lambda = 2, dims = 20, resolution = 0.1
  ),
  SPF_GF_FMT_Jej = list(
    treatment = c("SPF", "GF", "FMT"), segment = "Jej",
    lambda = 2, dims = 30, resolution = 0.1
  ),
  SPF_GF_FMT_Col = list(
    treatment = c("SPF", "GF", "FMT"), segment = "Col",
    lambda = 2, dims = 30, resolution = 0.1
  )
)

# cell type composition (LTSR), run on the SPF + GF + FMT objects
LTSR_SEGMENTS  <- c("Jej", "Col")
LTSR_GROUPS    <- c("SPF", "GF", "FMT")
LTSR_REFERENCE <- "GF"
LTSR_MAXFC     <- 3

# neighbourhood differential abundance (miloR)
# each comparison is c(reference, test); logFC > 0 means enriched in the test group
MILO_SEGMENTS    <- c("Jej", "Col")
MILO_COMPARISONS <- list(SPF_GF = c("GF", "SPF"), FMT_GF = c("GF", "FMT"))
MILO_REDUCTION   <- "harmony"
MILO_K           <- 30
MILO_D           <- 30
MILO_PROP        <- 0.1
MILO_FDR         <- 0.2

# differential expression between treatments, per cell type
DEG_SEGMENTS    <- c("Jej", "Col")
DEG_COMPARISONS <- list(SPF_GF = c("SPF", "GF"), FMT_GF = c("FMT", "GF"))
DEG_MIN_CELLS   <- 10
DEG_PADJ        <- 0.05
DEG_PCT         <- 0.3
DEG_LOG2FC      <- log2(1.5)

# cell-cell communication (CellChat)
# pig genes are renamed to human symbols because CellChatDB is human
ORTHOLOG_FILE <- file.path(PROJECT_DIR, "Sample_metdata", "pig2huam2mouse_ensembl106.csv")

CELLCHAT_SEGMENTS    <- c("Jej", "Col")
CELLCHAT_TREATMENTS  <- c("SPF", "GF", "FMT")
CELLCHAT_DB          <- "Secreted Signaling"
CELLCHAT_MIN_CELLS   <- 3
CELLCHAT_WORKERS     <- 4
CELLCHAT_COMPARISONS <- list(SPF_GF = c("SPF", "GF"), FMT_GF = c("FMT", "GF"))

# MetaNeighbor: reproducibility of the cell type labels
MN_SEGMENTS <- c("Jej", "Col")

# SCENIC: regulon activity, pySCENIC runs outside R between the two scripts
SCENIC_DIR     <- file.path(OUT_DIR, "scenic")
SCENIC_LOOM    <- file.path(SCENIC_DIR, "scenic.loom")
SCENIC_TF_LIST <- file.path(PROJECT_DIR, "Sample_metdata", "hs_hgnc_curated_tfs2pig.txt")
SCENIC_WORKERS <- 10
SCENIC_TOP     <- 10
