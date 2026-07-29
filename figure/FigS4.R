library(harmony)
library(Seurat)
library(ggplot2)
library(dplyr)
library(tibble)
library(data.table)
library(RColorBrewer)
library(Matrix)
library(patchwork)
library(Scillus)
library(tidyr)
library(ggpubr)

rm(list = ls()); gc()
ORIGINAL_DIR <- "D:/大学的资料/R/ESCC/data/FigS4"
setwd(ORIGINAL_DIR)


#--------S4A--------
tam <- readRDS(file = "D:/大学的资料/R/ESCC/data/TCN2_high_low/tam_group.RDS")
p1 <- FeaturePlot(
  tam,
  features = "TCN2",
  reduction = "umap",
  order = TRUE
) +
  labs(x = "UMAP1", y = "UMAP2",colour = "TCN2") +  
  ggtitle(NULL) +
  theme_classic() +
  theme(
    plot.title = element_text(size = 16, hjust = 0.5),
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 16),
    legend.text = element_text(size = 15),
    legend.title = element_text(size = 15)
  )
p1
ggsave(file.path("Fig7A-TCN2expression.png"), p1, width = 4.5, height = 4, dpi = 300, bg = "white")


#--------S4B--------
tam <- RenameIdents(
  tam,
  "TCN2-low" = "TCN2-low TAM_LYVE1",
  "TCN2-high" = "TCN2-high TAM_LYVE1"
)

meta <- tam@meta.data
meta$group <- Idents(tam)
plot_df <- meta %>%
  group_by(Sample, group) %>%
  summarise(
    n = n(),
    .groups = "drop"
  ) %>%
  group_by(Sample) %>%
  mutate(
    percent = n / sum(n) * 100
  )
p <- ggplot(
  plot_df,
  aes(
    x = Sample,
    y = percent,
    fill = group
  )
) +
  geom_bar(
    stat = "identity",
    width = 0.8
  ) +
  scale_fill_manual(
    name = "TCN2 group",
    values = c(
      "TCN2-low TAM_LYVE1" = "#A6CEE3",
      "TCN2-high TAM_LYVE1" = "#8B5A9E"
    )
  ) +
  scale_y_continuous(
    expand = c(0,0),
    labels = function(x) paste0(x, "%")
  ) +
  labs(
    x = NULL,
    y = "Cell proportion"
  ) +
  theme_classic() +
  theme(
    axis.text.x = element_text(
      angle = 45,
      hjust = 1,
      size = 11
    ),
    axis.text.y = element_text(
      size = 12
    ),
    axis.title.y = element_text(
      size = 14
    ),
    legend.title = element_text(
      size = 14,
      face = "bold"
    ),
    legend.text = element_text(
      size = 14
    ),
    
    legend.key.size = unit(
      0.6,
      "cm"
    )
  )

p
ggsave(file.path("FigS4B-TCN2expression.png"), p, width = 8, height = 4, dpi = 300, bg = "white")


#--------S4C--------
tam <- readRDS(file = "D:/大学的资料/R/ESCC/data/TCN2_high_low/tam_group_pos.RDS")
tam_1 <- readRDS(file = "D:/大学的资料/R/ESCC/data/TCN2_high_low/tam_group.RDS")

group_pos <- tam$tam_group
group_cluster <- tam_1$tam_group

same_cells <- all(
  names(group_pos) == names(group_cluster)
)

print(
  paste0(
    "Cell order identical: ",
    same_cells
  )
)

if(!same_cells){
  
  common_cells <- intersect(
    names(group_pos),
    names(group_cluster)
  )
  
  group_pos <- group_pos[common_cells]
  group_cluster <- group_cluster[common_cells]
}
conf_mat <- table(
  TCN2_positive_method = group_pos,
  Cluster_method = group_cluster
)

print(conf_mat)

accuracy <- sum(diag(conf_mat)) /
  sum(conf_mat)

print(
  paste0(
    "Accuracy = ",
    round(accuracy * 100, 2),
    "%"
  )
)

conf_prop_row <- prop.table(
  conf_mat,
  margin = 1
)

print(round(conf_prop_row, 3))

conf_prop_col <- prop.table(
  conf_mat,
  margin = 2
)

print(round(conf_prop_col, 3))

pdf(file = "FigS4C-heatmap.pdf", width = 4.3, height = 4)
p <- pheatmap(
  conf_mat,
  
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  
  display_numbers = TRUE,
  number_format = "%.0f",
  
  color = colorRampPalette(
    c(
      "#F7FBFF",
      "#C6DBEF",
      "#6BAED6",
      "#2171B5"
    )
  )(100),
  fontsize = 14,
  fontsize_number = 14,
  fontsize_row = 14,
  fontsize_col = 14,
  border_color = "grey80",
  main = NULL,
  legend = FALSE,
  legend_breaks = pretty(
    range(conf_mat),
    n = 4
  ),
  legend_labels = pretty(
    range(conf_mat),
    n = 4
  )
)
p
dev.off()


#--------S4D--------
tam <- readRDS(file = "D:/大学的资料/R/ESCC/data/TCN2_high_low/tam_group_pos.RDS")
options(stringsAsFactors = FALSE)
set.seed(20260713)

input_gmt <- "D:/大学的资料/R/ESCC/data/Fig7-new/c2.cp.reactome.v2026.1.Hs.symbols.gmt"

tam_group <- data.frame(
  barcode = colnames(tam),
  tam_group = tam$tam_group
)
tam_group <- tam_group[colnames(tam),]
tam$TCN2_group <- tam_group$tam_group
tam$TCN2_group <- factor(tam$TCN2_group, levels = c("TCN2-negative", "TCN2-positive"))
group_counts <- as.data.frame(table(tam$TCN2_group))
colnames(group_counts) <- c("group", "n_cells")

expr <- GetAssayData(tam, assay = "RNA", layer = "data")
de <- presto::wilcoxauc(expr, y = tam$TCN2_group)
de_positive <- de[de$group == "TCN2-positive", , drop = FALSE]
de_positive$avg_log2FC <- log2((de_positive$avgExpr + 1e-9) /
                             (de[de$group == "TCN2-negative", "avgExpr"] + 1e-9))

rank_keep <- de_positive$pct_in >= 1 | de_positive$pct_out >= 1
rank_stats <- 2 * (de_positive$auc[rank_keep] - 0.5) +
  de_positive$logFC[rank_keep] * 1e-7

names(rank_stats) <- de_positive$feature[rank_keep]
rank_stats <- sort(rank_stats[is.finite(rank_stats)], decreasing = TRUE)

pathways <- fgsea::gmtPathways(input_gmt)
gsea <- fgsea::fgseaMultilevel(
  pathways = pathways,
  stats = rank_stats,
  minSize = 10,
  maxSize = 500,
  eps = 0,
  nPermSimple = 10000
)
gsea <- as.data.frame(gsea)
if (anyNA(gsea$pval)) {
  na_paths <- gsea$pathway[is.na(gsea$pval)]
  fallback <- suppressWarnings(fgsea::fgseaSimple(
    pathways = pathways[na_paths], stats = rank_stats,
    nperm = 10000, minSize = 10, maxSize = 500
  ))
  fallback <- as.data.frame(fallback)
  replace_cols <- intersect(c("pval", "padj", "ES", "NES", "size", "leadingEdge"),
                            colnames(fallback))
  for (p in na_paths) {
    i <- match(p, gsea$pathway)
    j <- match(p, fallback$pathway)
    if (!is.na(j)) gsea[i, replace_cols] <- fallback[j, replace_cols]
  }
  gsea$padj <- p.adjust(gsea$pval, method = "BH")
}
gsea$enriched_in <- "TCN2-positive"
gsea$leadingEdge <- vapply(gsea$leadingEdge, paste, collapse = ";", FUN.VALUE = character(1))
gsea_positive <- gsea[order(gsea$padj, -abs(gsea$NES)), ]
gsea_sig <- gsea[!is.na(gsea$padj) & gsea$padj < 0.05, ]


#计算negative组
options(stringsAsFactors = FALSE)
set.seed(20260713)

input_gmt <- "D:/大学的资料/R/ESCC/data/Fig7-new/c2.cp.reactome.v2026.1.Hs.symbols.gmt"

tam_group <- data.frame(
  barcode = colnames(tam),
  tam_group = tam$tam_group
)
tam_group <- tam_group[colnames(tam),]
tam$TCN2_group <- tam_group$tam_group
tam$TCN2_group <- factor(tam$TCN2_group, levels = c("TCN2-negative", "TCN2-positive"))
group_counts <- as.data.frame(table(tam$TCN2_group))
colnames(group_counts) <- c("group", "n_cells")

expr <- GetAssayData(tam, assay = "RNA", layer = "data")
de <- presto::wilcoxauc(expr, y = tam$TCN2_group)
de_negative <- de[de$group == "TCN2-negative", , drop = FALSE]
de_negative$avg_log2FC <- log2((de_negative$avgExpr + 1e-9) /
                            (de[de$group == "TCN2-negative", "avgExpr"] + 1e-9))

rank_keep <- de_negative$pct_in >= 1 | de_negative$pct_out >= 1
rank_stats <- 2 * (de_negative$auc[rank_keep] - 0.5) +
  de_negative$logFC[rank_keep] * 1e-7

names(rank_stats) <- de_negative$feature[rank_keep]
rank_stats <- sort(rank_stats[is.finite(rank_stats)], decreasing = TRUE)

pathways <- fgsea::gmtPathways(input_gmt)
gsea <- fgsea::fgseaMultilevel(
  pathways = pathways,
  stats = rank_stats,
  minSize = 10,
  maxSize = 500,
  eps = 0,
  nPermSimple = 10000
)
gsea <- as.data.frame(gsea)
if (anyNA(gsea$pval)) {
  na_paths <- gsea$pathway[is.na(gsea$pval)]
  fallback <- suppressWarnings(fgsea::fgseaSimple(
    pathways = pathways[na_paths], stats = rank_stats,
    nperm = 10000, minSize = 10, maxSize = 500
  ))
  fallback <- as.data.frame(fallback)
  replace_cols <- intersect(c("pval", "padj", "ES", "NES", "size", "leadingEdge"),
                            colnames(fallback))
  for (p in na_paths) {
    i <- match(p, gsea$pathway)
    j <- match(p, fallback$pathway)
    if (!is.na(j)) gsea[i, replace_cols] <- fallback[j, replace_cols]
  }
  gsea$padj <- p.adjust(gsea$pval, method = "BH")
}
gsea$enriched_in <- "TCN2-negative"
gsea$leadingEdge <- vapply(gsea$leadingEdge, paste, collapse = ";", FUN.VALUE = character(1))
gsea_negative <- gsea[order(gsea$padj, -abs(gsea$NES)), ]
gsea_sig_negative <- gsea[!is.na(gsea$padj) & gsea$padj < 0.05, ]

pathway_keep <- c(
  "REACTOME_INTERFERON_GAMMA_SIGNALING",
  "REACTOME_INTERFERON_ALPHA_BETA_SIGNALING",
  "REACTOME_ANTIGEN_PROCESSING_CROSS_PRESENTATION",
  "REACTOME_MHC_CLASS_II_ANTIGEN_PRESENTATION",
  "REACTOME_ANTIGEN_PRESENTATION_FOLDING_ASSEMBLY_AND_PEPTIDE_LOADING_OF_CLASS_I_MHC",
  "REACTOME_IMMUNOREGULATORY_INTERACTIONS_BETWEEN_A_LYMPHOID_AND_A_NON_LYMPHOID_CELL",
  "REACTOME_INTERLEUKIN_10_SIGNALING"
)


pathway_order <- c(
  "REACTOME_INTERFERON_GAMMA_SIGNALING",
  "REACTOME_INTERFERON_ALPHA_BETA_SIGNALING",
  "REACTOME_ANTIGEN_PROCESSING_CROSS_PRESENTATION",
  "REACTOME_MHC_CLASS_II_ANTIGEN_PRESENTATION",
  "REACTOME_ANTIGEN_PRESENTATION_FOLDING_ASSEMBLY_AND_PEPTIDE_LOADING_OF_CLASS_I_MHC",
  "REACTOME_IMMUNOREGULATORY_INTERACTIONS_BETWEEN_A_LYMPHOID_AND_A_NON_LYMPHOID_CELL",
  "REACTOME_INTERLEUKIN_10_SIGNALING"
)

pathway_labels <- c(
  "REACTOME_INTERFERON_GAMMA_SIGNALING" =
    "IFN-gamma signaling",
  
  "REACTOME_INTERFERON_ALPHA_BETA_SIGNALING" =
    "IFN-alpha/beta signaling",
  
  "REACTOME_ANTIGEN_PROCESSING_CROSS_PRESENTATION" =
    "Cross-presentation",
  
  "REACTOME_MHC_CLASS_II_ANTIGEN_PRESENTATION" =
    "MHC-II antigen presentation",
  
  "REACTOME_ANTIGEN_PRESENTATION_FOLDING_ASSEMBLY_AND_PEPTIDE_LOADING_OF_CLASS_I_MHC" =
    "MHC-I antigen presentation",
  
  "REACTOME_IMMUNOREGULATORY_INTERACTIONS_BETWEEN_A_LYMPHOID_AND_A_NON_LYMPHOID_CELL" =
    "Lymphoid-myeloid immunoregulation",
  
  "REACTOME_INTERLEUKIN_10_SIGNALING" =
    "IL-10 signaling"
)
gsea <- rbind(gsea_positive, gsea_negative)
plot_df <- gsea %>%
  filter(pathway %in% pathway_keep) %>%
  dplyr::select(
    pathway,
    enriched_in,
    padj,
    NES
  ) %>%
  mutate(
    log10P = -log10(pmax(padj, 1e-300)),
    pathway = factor(
      pathway,
      levels = pathway_order
    ),
    enriched_in = factor(
      enriched_in,
      levels = c("TCN2-negative", "TCN2-positive")
    )
  )

# 检查筛选后的数据
print(plot_df)

p1 <- ggplot(
  plot_df,
  aes(
    x = pathway,
    y = enriched_in
  )
) +
  geom_point(
    aes(
      size = log10P,
      fill = NES
    ),
    shape = 21,
    color = "black",
    stroke = 0.6
  ) +
  scale_fill_gradient2(
    name = "NES",
    low = "#2166AC",
    mid = "white",
    high = "#B2182B",
    midpoint = 0,
    breaks = pretty_breaks(5)
  ) +
  scale_size_continuous(
    name = expression(-log[10](P[adj])),
    range = c(3, 10),
    breaks = pretty_breaks(4)
  ) +
  scale_x_discrete(
    position = "top",
    labels = pathway_labels,
    drop = FALSE
  ) +
  scale_y_discrete(
    position = "right",
    drop = FALSE
  ) +
  coord_cartesian(clip = "off") +
  theme_bw() +
  labs(
    x = NULL,
    y = NULL
  ) +
  theme(
    panel.border = element_rect(
      fill = NA,
      color = "black",
      linewidth = 1
    ),
    panel.grid.major = element_line(
      color = "grey90",
      linewidth = 0.4
    ),
    panel.grid.minor = element_blank(),
    
    axis.text.x = element_text(
      size = 12,
      angle = 55,
      colour = "black",
      hjust = 0,
      vjust = 0
    ),
    axis.text.y = element_text(
      size = 13,
      colour = "black"
    ),
    
    legend.title = element_text(size = 13),
    legend.text = element_text(size = 11),
    legend.position = "bottom",
    legend.box = "horizontal",
    
    plot.margin = margin(
      t = 5,
      r = 5,
      b = 2,
      l = 5,
      unit = "pt"
    )
  ) +
  guides(
    fill = guide_colorbar(
      title.position = "top",
      title.hjust = 0.5,
      order = 1,
      barwidth = unit(4, "cm")
    ),
    size = guide_legend(
      title.position = "top",
      title.hjust = 0.5,
      order = 2
    )
  )

p1
ggsave(file.path("FigS4D-heatmap.pdf"), p1, width =7, height = 4.3, bg = "white",
       device = cairo_pdf)
