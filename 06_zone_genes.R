# zone genes, one step: NMF programs -> zone genes (markers + SPARK) ->
# intersection -> hierarchical clustering, for the two SPF jejunum samples
source("00_config.R")
suppressPackageStartupMessages({
  library(Seurat)
  library(SPARK)
  library(NMF)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(ggplot2)
})

out_dir <- file.path(OUT_DIR, "zone")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

zone_pos <- function(z) as.numeric(sub("V", "", sub("Crypt", "V0", z))) + 1

read_zone <- function(id) {
  seu <- readRDS(file.path(OUT_DIR, "clean", paste0(id, ".rds")))
  seu <- subset(seu, cells = colnames(seu)[!is.na(seu$zone)])
  NormalizeData(seu, verbose = FALSE)
}

# 1. NMF programs, ordered along the crypt-villus axis
nmf_program <- function(seu, id) {
  dat <- FindVariableFeatures(seu, nfeatures = NMF_HVG, verbose = FALSE)
  dat <- ScaleData(dat, do.center = FALSE, verbose = FALSE)

  vm <- GetAssayData(dat, layer = "scale.data")
  vm <- vm[rowSums(vm) != 0, ]
  res <- nmf(vm, NMF_RANK)
  saveRDS(res, file.path(out_dir, paste0(id, "_nmf.rds")))

  coef <- t(NMF::coef(res))
  pos  <- zone_pos(seu$zone)
  coef <- coef[, order(apply(coef, 2, weighted.mean, w = pos))]
  colnames(coef) <- paste0("NMF_", seq_len(ncol(coef)))

  agg <- rowsum(coef, seu$zone)
  agg <- agg / as.vector(table(seu$zone)[rownames(agg)])

  heat <- as.data.frame(agg) %>%
    rownames_to_column("zone") %>%
    pivot_longer(-zone, names_to = "factor", values_to = "weight") %>%
    group_by(factor) %>%
    mutate(z = as.numeric(scale(weight))) %>%
    ungroup()

  lv <- names(sort(tapply(pos, seu$zone, mean)))

  p <- ggplot(heat, aes(factor, factor(zone, levels = lv), fill = z)) +
    geom_tile() +
    scale_fill_viridis_c(option = "B") +
    labs(x = NULL, y = "zone", title = id) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5))

  ggsave(file.path(out_dir, paste0(id, "_nmf_heatmap.pdf")), p, width = 6, height = 4)
}

# 2. zone genes by two methods
zone_markers <- function(seu, id) {
  Idents(seu) <- seu$zone
  markers <- FindAllMarkers(seu, only.pos = TRUE, min.pct = ZONE_MIN_PCT, verbose = FALSE)
  write.csv(markers, file.path(out_dir, paste0(id, "_deg.csv")), row.names = FALSE)

  counts <- GetAssayData(seu, layer = "counts")
  xy     <- GetTissueCoordinates(seu, cols = c("imagerow", "imagecol"))[colnames(seu), ]
  loc    <- as.matrix(xy[, 1:2])
  rownames(loc) <- colnames(counts)

  svg <- sparkx(counts, loc, numCores = SPARK_CORES, option = "mixture")
  write.csv(svg$res_mtest %>% rownames_to_column("gene"),
            file.path(out_dir, paste0(id, "_spark.csv")), row.names = FALSE)
}

# 3. keep genes that are both zone markers and spatially variable
select_genes <- function(id) {
  deg   <- read.csv(file.path(out_dir, paste0(id, "_deg.csv")))
  spark <- read.csv(file.path(out_dir, paste0(id, "_spark.csv")))

  keep <- intersect(
    deg$gene[deg$p_val_adj < ZONE_PADJ & deg$pct.1 > ZONE_MIN_PCT & deg$avg_log2FC > ZONE_LOG2FC],
    spark$gene[spark$adjustedPval < SPARK_PADJ]
  )
  keep <- sort(keep)

  write.csv(data.frame(gene = keep),
            file.path(out_dir, paste0(id, "_zone_genes.csv")), row.names = FALSE)
  keep
}

# 4. cluster the genes into co-varying programs
cluster_zone <- function(seu, id, genes) {
  bulk <- AverageExpression(seu, group.by = "zone", features = genes, layer = "data")[[1]]
  bulk <- t(apply(bulk, 1, function(x) (x - min(x)) / (max(x) - min(x))))

  cl  <- cutree(hclust(dist(bulk)), h = HCLUST_HEIGHT)
  tbl <- data.frame(gene = names(cl), cluster = as.integer(cl))
  write.csv(tbl, file.path(out_dir, paste0(id, "_zone_clusters.csv")), row.names = FALSE)

  df <- as.data.frame(bulk) %>%
    rownames_to_column("gene") %>%
    pivot_longer(-gene, names_to = "zone", values_to = "value") %>%
    left_join(tbl, by = "gene")
  df$zone <- factor(df$zone, levels = names(sort(tapply(zone_pos(df$zone), df$zone, mean))))

  means <- df %>% group_by(cluster, zone) %>% summarise(value = mean(value), .groups = "drop")

  p <- ggplot(df, aes(zone, value, group = gene)) +
    geom_line(colour = "grey80") +
    geom_line(data = means, aes(group = cluster), colour = "firebrick", linewidth = 1) +
    facet_wrap(~cluster, ncol = 4) +
    labs(x = "zone", y = "scaled expression", title = id) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5))

  ggsave(file.path(out_dir, paste0(id, "_zone_clusters.pdf")), p,
         height = length(unique(tbl$cluster)) * 0.6 + 2, width = 12)
}

for (id in ZONE_SAMPLES) {
  seu <- read_zone(id)
  nmf_program(seu, id)
  zone_markers(seu, id)
  cluster_zone(seu, id, select_genes(id))
}
