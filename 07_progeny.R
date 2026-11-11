# PROGENy pathway activity per spot
source("00_config.R")
suppressPackageStartupMessages({
  library(Seurat)
  library(progeny)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
})

out_dir <- file.path(OUT_DIR, "progeny")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

homo <- read.csv(ORTHOLOG_FILE, check.names = FALSE, stringsAsFactors = FALSE)
homo <- homo[homo[["Human homology type"]] == "ortholog_one2one" & !is.na(homo[["Human gene name"]]), ]
human_of <- tapply(homo[["Human gene name"]], homo[["Gene name"]], function(x) x[1])

progeny_one <- function(id) {
  seu <- readRDS(file.path(OUT_DIR, "clean", paste0(id, ".rds")))

  # pig symbols -> human symbols, one gene per human symbol
  mat <- GetAssayData(seu, layer = "counts")
  mat <- mat[!is.na(human_of[rownames(mat)]), ]
  mat <- mat[order(rowSums(mat), decreasing = TRUE), ]
  rownames(mat) <- human_of[rownames(mat)]
  mat <- mat[!duplicated(rownames(mat)), ]

  seu <- CreateSeuratObject(counts = mat, meta.data = seu@meta.data, assay = "Spatial")
  seu <- NormalizeData(seu, verbose = FALSE)

  act <- progeny(seu, scale = FALSE, assay_name = "Spatial", organism = "Human",
                 top = PROGENY_TOP, return_assay = TRUE, z_scores = FALSE)
  DefaultAssay(act) <- "progeny"
  act <- ScaleData(act, verbose = FALSE)

  scores <- as.data.frame(t(GetAssayData(act, layer = "scale.data", assay = "progeny")))
  scores$cell      <- rownames(scores)
  scores$sample    <- id
  scores$treatment <- unique(seu$treatment)
  scores$segment   <- unique(seu$segment)
  scores$zone      <- seu$zone[match(scores$cell, colnames(seu))]
  scores
}

scores <- bind_rows(lapply(unique(read.csv(METADATA_FILE)$Sample), progeny_one))
write.csv(scores, file.path(out_dir, "scores.csv"), row.names = FALSE)

means <- scores %>%
  pivot_longer(cols = -c(cell, sample, treatment, segment, zone),
               names_to = "pathway", values_to = "score") %>%
  group_by(treatment, segment, pathway) %>%
  summarise(score = mean(score), .groups = "drop")

p <- ggplot(means, aes(segment, pathway, fill = score)) +
  geom_tile() +
  facet_wrap(~treatment) +
  scale_fill_viridis_c() +
  labs(x = NULL, y = NULL) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5))

ggsave(file.path(out_dir, "progeny_heatmap.pdf"), p, width = 8, height = 5)
