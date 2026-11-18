# miloR: neighbourhood differential abundance between treatments, per segment
source("00_config.R")
suppressPackageStartupMessages({
  library(Seurat)
  library(miloR)
  library(SingleCellExperiment)
  library(dplyr)
  library(ggplot2)
  library(ggforce)
})

out_dir <- file.path(OUT_DIR, "milo")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

milo_one <- function(segment, comparison) {
  groups <- MILO_COMPARISONS[[comparison]]
  ref <- groups[1]
  test <- groups[2]

  seu <- readRDS(file.path(OUT_DIR, "merge", paste0("SPF_GF_FMT_", segment, ".rds")))
  seu <- subset(seu, cells = colnames(seu)[seu$treatment %in% groups])
  seu$treatment <- factor(seu$treatment, levels = groups)

  milo <- Milo(as.SingleCellExperiment(seu))
  milo <- buildGraph(milo, k = MILO_K, d = MILO_D, reduced.dim = MILO_REDUCTION)
  milo <- makeNhoods(milo, prop = MILO_PROP, k = MILO_K, d = MILO_D, refined = TRUE)
  milo <- countCells(milo, meta.data = as.data.frame(colData(milo)), sample = "sample")
  milo <- calcNhoodDistance(milo, d = MILO_D, reduced.dim = MILO_REDUCTION)

  design <- colData(milo) %>%
    as.data.frame() %>%
    distinct(sample, .keep_all = TRUE) %>%
    select(sample, treatment)
  rownames(design) <- colnames(milo@nhoodCounts)

  da <- testNhoods(milo, design = ~ treatment, design.df = design)
  milo <- buildNhoodGraph(milo)

  da <- annotateNhoods(milo, da, coldata_col = "cell_type")
  da$group <- ifelse(da$SpatialFDR < MILO_FDR & da$logFC > 0, test,
                     ifelse(da$SpatialFDR < MILO_FDR & da$logFC < 0, ref, "ns"))
  da$group <- factor(da$group, levels = c(test, ref, "ns"))

  write.csv(da, file.path(out_dir, paste0(segment, "_", comparison, "_da.csv")), row.names = FALSE)

  p_graph <- plotNhoodGraphDA(milo, da, alpha = 0.1) +
    scale_fill_gradientn(colours = c("#2166AC", "#F7F7F7", "#B2182B"))

  p_sina <- ggplot(da, aes(logFC, cell_type, colour = group)) +
    geom_sina(alpha = 0.5, method = "density", scale = "width", maxwidth = 0.8) +
    scale_colour_manual(values = c("#B2182B", "#2166AC", "grey70")) +
    labs(x = paste0("logFC (", test, " vs ", ref, ")"), y = NULL,
         title = paste(segment, comparison)) +
    theme_bw()

  ggsave(file.path(out_dir, paste0(segment, "_", comparison, "_graph.pdf")), p_graph, width = 6, height = 5)
  ggsave(file.path(out_dir, paste0(segment, "_", comparison, "_sina.pdf")),  p_sina,  width = 5, height = 4)
}

for (segment in MILO_SEGMENTS) {
  for (comparison in names(MILO_COMPARISONS)) {
    milo_one(segment, comparison)
  }
}
