# regulon activity from the pySCENIC loom: AUC, RSS, heatmap
source("00_config.R")
suppressPackageStartupMessages({
  library(Seurat)
  library(SCopeLoomR)
  library(SCENIC)
  library(AUCell)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(pheatmap)
})

seu <- readRDS(file.path(SCENIC_DIR, "seurat_scenic.rds"))

loom <- open_loom(SCENIC_LOOM)
auc  <- get_regulons_AUC(loom, column.attr.name = "RegulonsAUC")
close_loom(loom)

# regulons x cells, cell annotation in the same order
auc <- getAUC(auc)
ann <- setNames(seu$group_label, colnames(seu))[colnames(auc)]

rss <- calcRSS(AUC = auc, cellAnnotation = ann)
write.csv(rss, file.path(SCENIC_DIR, "rss.csv"))

# top regulons per group
long <- as.data.frame(rss) %>%
  tibble::rownames_to_column("regulon") %>%
  pivot_longer(-regulon, names_to = "group", values_to = "rss")

top <- long %>%
  group_by(group) %>%
  slice_max(rss, n = SCENIC_TOP) %>%
  ungroup()

p <- ggplot(top, aes(rss, reorder(regulon, rss))) +
  geom_col(fill = "grey40") +
  facet_wrap(~group, scales = "free_y", ncol = 4) +
  labs(x = "RSS", y = NULL) +
  theme_bw() +
  theme(axis.text.y = element_text(size = 6))

ggsave(file.path(SCENIC_DIR, "rss_top.pdf"), p, width = 14, height = 10)

# regulon-regulon correlation
pdf(file.path(SCENIC_DIR, "rss_heatmap.pdf"), width = 6, height = 6)
pheatmap(cor(rss), border_color = NA, treeheight_row = 60, treeheight_col = 60, fontsize = 4)
dev.off()
