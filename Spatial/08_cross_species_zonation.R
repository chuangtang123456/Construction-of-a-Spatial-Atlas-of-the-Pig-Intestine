# cross-species zonation: pig vs human vs mouse expression along the villus axis
source("00_config.R")
suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(ggplot2)
})

out_dir <- file.path(OUT_DIR, "cross_species")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# pig -> human symbols, one-to-one in both species
homo <- read.csv(ORTHOLOG_FILE, check.names = FALSE, stringsAsFactors = FALSE)
homo <- homo[homo[["Human homology type"]] == "ortholog_one2one" &
               homo[["Mouse homology type"]] == "ortholog_one2one" &
               !is.na(homo[["Human gene name"]]), ]
human_of <- tapply(homo[["Human gene name"]], homo[["Gene name"]], function(x) x[1])

# pig: mean expression per zone, averaged over samples
pig_one <- function(id) {
  seu <- readRDS(file.path(OUT_DIR, "clean", paste0(id, ".rds")))
  seu <- subset(seu, cells = colnames(seu)[!is.na(seu$zone)])
  seu <- NormalizeData(seu, verbose = FALSE)

  bulk <- AverageExpression(seu, group.by = "zone", layer = "data")[[1]]
  colnames(bulk) <- sub("^Spatial\\.", "", colnames(bulk))
  bulk[, grepl("^V[0-9]+$", colnames(bulk)), drop = FALSE]
}

mats  <- lapply(CROSS_SAMPLES, pig_one)
genes <- Reduce(intersect, lapply(mats, rownames))
pig   <- Reduce(`+`, lapply(mats, function(m) m[genes, ])) / length(mats)

pig <- pig[!is.na(human_of[rownames(pig)]), ]
pig <- pig[order(rowSums(pig), decreasing = TRUE), ]
rownames(pig) <- human_of[rownames(pig)]
pig <- pig[!duplicated(rownames(pig)), ]
colnames(pig) <- paste("Pig", "mean", colnames(pig), sep = "_")

# human / mouse reference
hm <- read.csv(CROSS_SPECIES_REF, skip = 1, row.names = 1)
hm <- hm[, grep("mean", colnames(hm))]
hm <- hm[apply(hm[, grep("Human", colnames(hm))], 1, max) > CROSS_MIN_EXPR &
           apply(hm[, grep("Mouse", colnames(hm))], 1, max) > CROSS_MIN_EXPR, ]

keep <- intersect(rownames(pig), rownames(hm))
all  <- cbind(pig[keep, ], hm[keep, ])

# expression scaled to [0, 1] per gene and species, zones scaled to [0, 1]
long <- all %>%
  rownames_to_column("gene") %>%
  pivot_longer(-gene) %>%
  mutate(species = sub("_.*$", "", name),
         zone    = as.numeric(sub("^V", "", sub(".*_", "", name)))) %>%
  group_by(species) %>%
  mutate(x = (zone - min(zone)) / (max(zone) - min(zone))) %>%
  ungroup() %>%
  group_by(gene, species) %>%
  mutate(y = (value - min(value)) / (max(value) - min(value))) %>%
  ungroup()

# curve similarity per gene and species pair
grid <- seq(0, 1, length.out = 101)

curve_stats <- function(df) {
  cur <- df %>%
    arrange(species, x) %>%
    group_by(species) %>%
    summarise(y = list(approx(x, y, xout = grid, rule = 2, ties = mean)$y), .groups = "drop")
  if (nrow(cur) < 2) return(NULL)

  m  <- do.call(rbind, cur$y)
  sp <- cur$species
  out <- list()
  for (i in seq_len(nrow(cur) - 1)) {
    for (j in seq(i + 1, nrow(cur))) {
      out[[length(out) + 1]] <- tibble(
        pair     = paste(sp[i], sp[j], sep = "-"),
        pearson  = cor(m[i, ], m[j, ]),
        center_1 = sum(grid * m[i, ]) / sum(m[i, ]),
        center_2 = sum(grid * m[j, ]) / sum(m[j, ])
      )
    }
  }
  bind_rows(out)
}

res <- long %>%
  group_by(gene) %>%
  group_modify(~ curve_stats(.x)) %>%
  ungroup() %>%
  mutate(center_diff = center_1 - center_2,
         status = ifelse(abs(center_diff) > CROSS_CENTER_CUT & pearson < CROSS_PEARSON_CUT,
                         "Discordant", "Concordant"))

write.csv(res, file.path(out_dir, "curve_similarity.csv"), row.names = FALSE)

p <- ggplot(res, aes(center_diff, pearson, colour = status)) +
  geom_point(size = 0.6, alpha = 0.6) +
  facet_wrap(~pair, ncol = 1) +
  labs(x = "centre of mass shift", y = "pearson r") +
  theme_bw() +
  theme(legend.position = "top")

ggsave(file.path(out_dir, "curve_similarity.pdf"), p, width = 4, height = 7)
