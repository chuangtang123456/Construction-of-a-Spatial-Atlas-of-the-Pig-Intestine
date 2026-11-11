# spatial aggregation: how many same-type spots sit among the k nearest neighbours
source("00_config.R")
suppressPackageStartupMessages({
  library(Seurat)
  library(FNN)
  library(dplyr)
  library(ggplot2)
})

meta <- read.csv(METADATA_FILE, stringsAsFactors = FALSE)

aggregation_one <- function(id) {
  w   <- read.csv(file.path(OUT_DIR, "rctd", paste0(id, ".csv")), row.names = 1, check.names = FALSE)
  seu <- readRDS(file.path(OUT_DIR, "clean", paste0(id, ".rds")))

  xy <- GetTissueCoordinates(seu, cols = c("imagerow", "imagecol"))[rownames(w), ]
  nn <- get.knn(as.matrix(xy[, 1:2]), k = CELL_DIST_K + 1)$nn.index[, -1]

  rows <- lapply(colnames(w), function(ct) {
    signal <- w[[ct]]
    sel    <- signal > quantile(signal, CELL_DIST_QUANTILE, na.rm = TRUE)
    k      <- rowSums(matrix(sel[nn], nrow = nrow(nn)))[sel]
    if (length(k) < 2) return(NULL)

    data.frame(cell_type = ct, n_same = 0:CELL_DIST_K,
               n = as.integer(table(factor(k, levels = 0:CELL_DIST_K))))
  })

  dist <- bind_rows(rows) %>%
    group_by(cell_type) %>%
    mutate(pct = n / sum(n)) %>%
    ungroup()

  out_dir <- file.path(OUT_DIR, "cell_aggregation")
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  write.csv(dist, file.path(out_dir, paste0(id, ".csv")), row.names = FALSE)

  p <- ggplot(dist, aes(n_same, pct)) +
    geom_col(fill = "grey40", width = 0.8) +
    facet_wrap(~cell_type, ncol = 6, scales = "free_y") +
    labs(x = "same-type neighbours (of 6)", y = "proportion", title = id) +
    theme_bw() +
    theme(panel.grid.minor = element_blank())

  ggsave(file.path(out_dir, paste0(id, ".pdf")), p, width = 12, height = 6)
}

for (id in unique(meta$Sample)) aggregation_one(id)
