# differential expression between treatments, per cell type and segment
source("00_config.R")
suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
})

out_dir <- file.path(OUT_DIR, "deg")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

deg_one <- function(segment) {
  seu <- readRDS(file.path(OUT_DIR, "merge", paste0("SPF_GF_FMT_", segment, ".rds")))
  seu$for_deg <- paste0(seu$group, "_", seu$cell_type)
  Idents(seu) <- seu$for_deg

  # a cell type enters the test only if every treatment has enough cells
  n_by_group <- table(seu$cell_type, seu$treatment)
  cell_types <- rownames(n_by_group)[apply(n_by_group, 1, min) >= DEG_MIN_CELLS]

  res <- list()
  for (cmp in names(DEG_COMPARISONS)) {
    trt <- DEG_COMPARISONS[[cmp]]
    for (ct in cell_types) {
      d <- FindMarkers(seu,
                       ident.1 = paste(trt[1], segment, ct, sep = "_"),
                       ident.2 = paste(trt[2], segment, ct, sep = "_"),
                       verbose = FALSE)
      if (nrow(d) == 0) next

      d$gene       <- rownames(d)
      d$cell_type  <- ct
      d$comparison <- cmp
      res[[paste(cmp, ct)]] <- d
    }
  }

  degs <- bind_rows(res)
  write.csv(degs, file.path(out_dir, paste0(segment, "_all.csv")), row.names = FALSE)

  clean <- degs %>%
    filter(p_val_adj < DEG_PADJ,
           pct.1 > DEG_PCT | pct.2 > DEG_PCT,
           abs(avg_log2FC) > DEG_LOG2FC)

  # a gene is shared when it moves in the same direction in both comparisons
  wide <- clean %>%
    select(gene, cell_type, comparison, avg_log2FC) %>%
    pivot_wider(names_from = comparison, values_from = avg_log2FC) %>%
    mutate(shared = !is.na(SPF_GF) & !is.na(FMT_GF) & sign(SPF_GF) == sign(FMT_GF))

  clean <- clean %>%
    left_join(select(wide, gene, cell_type, shared), by = c("gene", "cell_type"))

  write.csv(clean, file.path(out_dir, paste0(segment, "_clean.csv")), row.names = FALSE)

  counts <- clean %>% count(cell_type, comparison)

  p <- ggplot(counts, aes(reorder(cell_type, -n), n, fill = comparison)) +
    geom_col(position = position_dodge2(), colour = "black") +
    scale_y_continuous(expand = c(0, 0)) +
    labs(x = NULL, y = "DEGs", title = segment) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))

  ggsave(file.path(out_dir, paste0(segment, "_counts.pdf")), p, width = 7, height = 4)
}

for (segment in DEG_SEGMENTS) deg_one(segment)
