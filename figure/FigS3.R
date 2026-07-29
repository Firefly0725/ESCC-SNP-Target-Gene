library(Seurat)
library(dplyr)
library(ggplot2)
library(dplyr)
library(tidyr)
library(ggpubr)

rm(list = ls()); gc()
ORIGINAL_DIR <- " "
setwd(ORIGINAL_DIR)
scdata <- readRDS(file = "scdata.RDS")

#--------FigA--------
p1 <- DimPlot(scdata, 
              group.by = "cluster", 
              label = TRUE, 
              label.size = 5) +
  labs(x = "UMAP1", y = "UMAP2") +  
  theme(
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 14),
    legend.text = element_text(size = 14),
    plot.title = element_blank()
  )+
  guides(color = guide_legend(
    ncol = 2,
    override.aes = list(size = 5)     
  ))
p1
ggsave(file.path("Fig3SA_cluster_UMAP.png"), p1, width = 5.5, height = 4, dpi = 300, bg = "white")

#--------FigB--------
p1 <- FeaturePlot(
  object = scdata,
  features = "TAM_LYVE1_score1",
  reduction = "umap",
  order = TRUE,
  min.cutoff = "q05",
  max.cutoff = "q95"
) +
  labs(x = "UMAP1", y = "UMAP2",colour = "TAM_LYVE1\nScore",
       title = NULL) +
  theme_classic() +
  theme(
    plot.title = element_text(size = 14, hjust = 0.5),
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 14),
    legend.text = element_text(size = 14),
    legend.title = element_text(size = 14)
  )
p1
ggsave(file.path("Fig3SB_cluster_score.png"), p1, width = 5.5, height = 4, dpi = 300, bg = "white")

feature_plot <- function(genename){
  p1 <- FeaturePlot(
    object = scdata,
    features = genename,
    reduction = "umap",
    order = TRUE,
    min.cutoff = "q05",
    max.cutoff = "q95"
  ) +
    labs(x = NULL, y = NULL, colour = genename, title = NULL) +  
    theme_classic() +
    theme(
      plot.title = element_text(size = 14, hjust = 0.5),
      axis.title = element_blank(),           
      axis.text = element_blank(),      
      axis.ticks = element_blank(),        
      axis.line = element_blank(),        
      legend.text = element_text(size = 14),
      legend.title = element_text(size = 14)
    )
  return(p1)
}
p1 <- feature_plot("C1QB")
p2 <- feature_plot("CCL18")
p3 <- feature_plot("TREM2")
p4 <- feature_plot("CD209")
combined_plot <- p1 | p2 | p3 | p4
combined_plot
ggsave(
  filename = "Fig3SB_foursamples.png", 
  plot = combined_plot,
  width = 12.5,
  height = 3,
  dpi = 300
)


#--------FigC--------
cluster_score <- scdata@meta.data %>%
  group_by(cluster) %>%
  summarise(
    mean_score = mean(TAM_LYVE1_score1, na.rm = TRUE),
    median_score = median(TAM_LYVE1_score1, na.rm = TRUE),
    n_cell = n(),
    .groups = "drop"
  ) %>%
  arrange(desc(mean_score))
p1 <- ggplot(
  cluster_score,
  aes(
    x = mean_score,
    y = reorder(as.character(cluster), mean_score, decreasing = T)
  )
) +
  geom_col(fill = "blue") + 
  coord_flip() +
  theme_classic() +
  labs(
    y = "Cluster",
    x = "Mean TAM_LYVE1 module score"
  ) +
  theme(
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 14)
  )
p1

p1 <- ggplot(
  cluster_score,
  aes(
    x = mean_score,
    y = reorder(as.character(cluster), mean_score, decreasing = T)
  )
) +
  geom_col(fill = "blue") + 
  coord_flip() +
  theme_classic() +
  labs(
    y = "Cluster",
    x = "Mean TAM_LYVE1 module score"
  ) +
  theme(
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 14)
  )
p1
ggsave(file.path("Fig3SC-1_barPlot.pdf"), p1, width = 6.5, height = 4, bg = "white")

ordered_clusters <- cluster_score$cluster[order(cluster_score$mean_score, decreasing = TRUE)]
p1 <- VlnPlot(
  object = scdata,
  features = "TAM_LYVE1_score1",
  group.by = "cluster",
  pt.size = 0
) +
  ggtitle(NULL) +
  xlab("Cluster") +
  ylab("Module score") +
  scale_x_discrete(limits = ordered_clusters) + 
  theme_classic() +
  theme(
    plot.title = element_text(size = 16, hjust = 0.5),
    axis.title = element_text(size = 16),
    axis.text.x = element_text(
      size = 14,
      angle = 45,
      hjust = 1
    ),
    axis.text.y = element_text(size = 14),
    legend.position = "none"
  )
p1
ggsave(file.path("Fig3SC-1_BoxPlot.pdf"), p1, width = 6.5, height = 4, bg = "white")
