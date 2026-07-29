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
library(ggplot2)
library(dplyr)
library(RColorBrewer)
library(ComplexHeatmap)
library(circlize)
library(stats)
library(grid)
library(Seurat)
library(tidyr)
library(arrow)
library(scales)
library(ggvenn)
library(fgsea)
library(data.table)
library(tibble)
library(cowplot)
library(ggpubr)
library(patchwork)
library(ggforce)

rm(list = ls()); gc()
ORIGINAL_DIR <- "Fig7-new"
setwd(ORIGINAL_DIR)
tam <- readRDS(file = "tam_group.RDS")

#--------7A--------
tam <- RenameIdents(
  tam,
  "TCN2-low" = "TCN2−low TAM_LYVE1",
  "TCN2-high" = "TCN2-high TAM_LYVE1"
)

my_colors <- c(
  "TCN2−low TAM_LYVE1" = "#A6CEE3",
  "TCN2-high TAM_LYVE1" = "#8B5A9E"
)

p1 <- DimPlot(
  tam,
  cols = my_colors,
  label.size = 5
) +
  labs(
    x = "UMAP1",
    y = "UMAP2"
  ) +
  theme(
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 14),
    legend.text = element_text(size = 14),
    plot.title = element_blank(),
    legend.position = "bottom",
    legend.direction = "horizontal"
  ) +
  guides(
    color = guide_legend(
      ncol = 1,
      override.aes = list(size = 5)
    )
  )

p1

ggsave(
  file.path("Fig7A-UMAP.png"),
  p1,
  width = 3.5,
  height = 4,
  dpi = 300,
  bg = "white"
)

#--------7B--------
options(stringsAsFactors = FALSE)
set.seed(20260713)

input_gmt <- "h.all.v2026.1.Hs.symbols.gmt"

tam_group <- data.frame(
  barcode = colnames(tam),
  tam_group = tam$tam_group
)
tam_group <- tam_group[colnames(tam),]
tam$TCN2_group <- tam_group$tam_group
tam$TCN2_group <- factor(tam$TCN2_group, levels = c("TCN2-low", "TCN2-high"))
group_counts <- as.data.frame(table(tam$TCN2_group))
colnames(group_counts) <- c("group", "n_cells")

expr <- GetAssayData(tam, assay = "RNA", layer = "data")
de <- presto::wilcoxauc(expr, y = tam$TCN2_group)
de_high <- de[de$group == "TCN2-high", , drop = FALSE]
de_high$avg_log2FC <- log2((de_high$avgExpr + 1e-9) /
                             (de[de$group == "TCN2-low", "avgExpr"] + 1e-9))

rank_keep <- de_high$pct_in >= 1 | de_high$pct_out >= 1
rank_stats <- 2 * (de_high$auc[rank_keep] - 0.5) +
  de_high$logFC[rank_keep] * 1e-7

names(rank_stats) <- de_high$feature[rank_keep]
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
gsea$enriched_in <- "TCN2-high"
gsea$leadingEdge <- vapply(gsea$leadingEdge, paste, collapse = ";", FUN.VALUE = character(1))
gsea_high <- gsea[order(gsea$padj, -abs(gsea$NES)), ]
gsea_sig <- gsea[!is.na(gsea$padj) & gsea$padj < 0.05, ]

options(stringsAsFactors = FALSE)
set.seed(20260713)

input_gmt <- "h.all.v2026.1.Hs.symbols.gmt"

tam_group <- data.frame(
  barcode = colnames(tam),
  tam_group = tam$tam_group
)
tam_group <- tam_group[colnames(tam),]
tam$TCN2_group <- tam_group$tam_group
tam$TCN2_group <- factor(tam$TCN2_group, levels = c("TCN2-low", "TCN2-high"))
group_counts <- as.data.frame(table(tam$TCN2_group))
colnames(group_counts) <- c("group", "n_cells")

expr <- GetAssayData(tam, assay = "RNA", layer = "data")
de <- presto::wilcoxauc(expr, y = tam$TCN2_group)
de_low <- de[de$group == "TCN2-low", , drop = FALSE]
de_low$avg_log2FC <- log2((de_low$avgExpr + 1e-9) /
                            (de[de$group == "TCN2-low", "avgExpr"] + 1e-9))

rank_keep <- de_low$pct_in >= 1 | de_low$pct_out >= 1
rank_stats <- 2 * (de_low$auc[rank_keep] - 0.5) +
  de_low$logFC[rank_keep] * 1e-7

names(rank_stats) <- de_low$feature[rank_keep]
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
gsea$enriched_in <- "TCN2-low"
gsea$leadingEdge <- vapply(gsea$leadingEdge, paste, collapse = ";", FUN.VALUE = character(1))
gsea_low <- gsea[order(gsea$padj, -abs(gsea$NES)), ]
gsea_sig_low <- gsea[!is.na(gsea$padj) & gsea$padj < 0.05, ]

write.csv(gsea_sig, file = "gsea_result.csv", row.names = F)

pathway_keep <- c(
  "HALLMARK_INTERFERON_GAMMA_RESPONSE",
  "HALLMARK_INTERFERON_ALPHA_RESPONSE",
  "HALLMARK_ALLOGRAFT_REJECTION",
  "HALLMARK_INFLAMMATORY_RESPONSE",
  "HALLMARK_IL6_JAK_STAT3_SIGNALING",
  "HALLMARK_COMPLEMENT",

  "HALLMARK_HYPOXIA",
  "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION",
  "HALLMARK_OXIDATIVE_PHOSPHORYLATION",
  "HALLMARK_ANGIOGENESIS"
)

pathway_order <- c(
  "HALLMARK_INTERFERON_GAMMA_RESPONSE",
  "HALLMARK_INTERFERON_ALPHA_RESPONSE",
  "HALLMARK_ALLOGRAFT_REJECTION",
  "HALLMARK_INFLAMMATORY_RESPONSE",
  "HALLMARK_IL6_JAK_STAT3_SIGNALING",
  "HALLMARK_COMPLEMENT",
  "HALLMARK_HYPOXIA",
  "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION",
  "HALLMARK_OXIDATIVE_PHOSPHORYLATION",
  "HALLMARK_ANGIOGENESIS"
)

pathway_labels <- c(
  "HALLMARK_INTERFERON_GAMMA_RESPONSE" =
    "IFN-\u03b3 response",
  "HALLMARK_INTERFERON_ALPHA_RESPONSE" =
    "IFN-\u03b1 response",
  "HALLMARK_ALLOGRAFT_REJECTION" =
    "Allograft rejection",
  "HALLMARK_INFLAMMATORY_RESPONSE" =
    "Inflammatory response",
  "HALLMARK_IL6_JAK_STAT3_SIGNALING" =
    "IL6–JAK–STAT3 signaling",
  "HALLMARK_COMPLEMENT" =
    "Complement",
  "HALLMARK_HYPOXIA" =
    "Hypoxia",
  "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION" =
    "Epithelial–mesenchymal transition",
  "HALLMARK_OXIDATIVE_PHOSPHORYLATION" =
    "Oxidative phosphorylation",
  "HALLMARK_ANGIOGENESIS" =
    "Angiogenesis"
)
gsea <- rbind(gsea_high, gsea_low)
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
      levels = c("TCN2-low", "TCN2-high")
    )
  )

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
ggsave(file.path("Fig7B-heatmap.pdf"), p1, width =8, height = 4.3, bg = "white",
       device = cairo_pdf)

#--------7C--------
tam <- readRDS(file = "tam_group.RDS")
options(stringsAsFactors = FALSE)
set.seed(20260713)

input_gmt <- "c2.cp.reactome.v2026.1.Hs.symbols.gmt"

tam_group <- data.frame(
  barcode = colnames(tam),
  tam_group = tam$tam_group
)
tam_group <- tam_group[colnames(tam),]
tam$TCN2_group <- tam_group$tam_group
tam$TCN2_group <- factor(tam$TCN2_group, levels = c("TCN2-low", "TCN2-high"))
group_counts <- as.data.frame(table(tam$TCN2_group))
colnames(group_counts) <- c("group", "n_cells")

expr <- GetAssayData(tam, assay = "RNA", layer = "data")
de <- presto::wilcoxauc(expr, y = tam$TCN2_group)
de_high <- de[de$group == "TCN2-high", , drop = FALSE]
de_high$avg_log2FC <- log2((de_high$avgExpr + 1e-9) /
                             (de[de$group == "TCN2-low", "avgExpr"] + 1e-9))

rank_keep <- de_high$pct_in >= 1 | de_high$pct_out >= 1
rank_stats <- 2 * (de_high$auc[rank_keep] - 0.5) +
  de_high$logFC[rank_keep] * 1e-7

names(rank_stats) <- de_high$feature[rank_keep]
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
gsea$enriched_in <- ifelse(gsea$ES > 0, "TCN2-high", "TCN2-low")
gsea$leadingEdge <- vapply(gsea$leadingEdge, paste, collapse = ";", FUN.VALUE = character(1))
gsea <- gsea[order(gsea$padj, -abs(gsea$NES)), ]
gsea_sig <- gsea[!is.na(gsea$padj) & gsea$padj < 0.05, ]
write.csv(gsea_sig, file = "gsea_result.csv",row.names = F)

plot_draw <- function(pathway_name, short_name, location){
  pathway_genes <- pathways[[pathway_name]]

  plotEnrichmentData <- function(pathway, stats,
                                 gseaParam=1) {

    if (any(!is.finite(stats))){
      stop("Not all stats values are finite numbers")
    }

    rnk <- rank(-stats)
    ord <- order(rnk)

    statsAdj <- stats[ord]
    statSigns <- sign(statsAdj)
    statsAdj <- (abs(statsAdj) ^ gseaParam)

    pathway <- unname(as.vector(na.omit(match(pathway, names(statsAdj)))))
    pathway <- sort(pathway)
    pathway <- unique(pathway)

    gseaRes <- calcGseaStat(statsAdj, selectedStats = pathway,
                            returnAllExtremes = TRUE)

    bottoms <- gseaRes$bottoms
    tops <- gseaRes$tops

    n <- length(statsAdj)
    xs <- as.vector(rbind(pathway - 1, pathway))
    ys <- as.vector(rbind(bottoms, tops))
    toPlot <- data.table(rank=c(0, xs, n + 1), ES=c(0, ys, 0))
    ticks <- data.table(rank=pathway, stat=statSigns[pathway]*statsAdj[pathway])
    stats <- data.table(rank=seq_along(stats), stat=statSigns*statsAdj)

    res <- list(
      curve=toPlot,
      ticks=ticks,
      stats=stats,
      posES=max(tops),
      negES=min(bottoms),
      spreadES=max(tops)-min(bottoms),
      maxAbsStat=max(abs(statsAdj)))
  }
  pd <- plotEnrichmentData(pathway_genes, rank_stats, gseaParam = 1)
  p <- ggplot(pd$curve, aes(x = rank, y = ES)) +
    geom_line(color = "seagreen", linewidth = 1) +
    geom_segment(data = pd$ticks,
                 aes(x = rank, y = -pd$spreadES/16,
                     xend = rank, yend = pd$spreadES/16),
                 linewidth = 0.2) +
    geom_hline(yintercept = 0, colour = "black") +
    labs(title = short_name,
         x = "Rank in ordered gene list",
         y = "Geneset enrichment score") +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5, size = 14),
          axis.text = element_text(color = "black", size = 13),
          axis.title = element_text(size = 14,color = "black"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          axis.line = element_line(color = "black", linewidth = 0.5),
          axis.ticks = element_line(color = "black"))

  res_line <- gsea_sig[gsea_sig$pathway == pathway_name, ]
  nes_val <- round(res_line$NES, 3)
  padj_val <- res_line$padj
  p_label <- ifelse(padj_val < 0.001, "<0.001", sprintf("%.3f", padj_val))
  p <- p + annotate("text",
                    x = Inf, y = location,
                    label = paste0("NES = ", nes_val, "\nPadj ", p_label),
                    hjust = 1.1, vjust = 0,
                    size = 5)
  p
  return(p)
}

p1 <- plot_draw(pathway_name = "REACTOME_INTERFERON_GAMMA_SIGNALING",
                short_name = "IFN-γ signaling",
                location = 0.6)
p1
p2 <- plot_draw(pathway_name = "REACTOME_MHC_CLASS_II_ANTIGEN_PRESENTATION",
                short_name = "MHC-II antigen presentation",
                location = 0.52)
p2
p3 <- plot_draw(pathway_name = "REACTOME_INTERLEUKIN_10_SIGNALING",
                short_name = "IL-10 signaling",
                location = 0.6)
p3
p4 <- plot_draw(pathway_name = "REACTOME_FCGR3A_MEDIATED_IL10_SYNTHESIS",
                short_name = "FCGR3A-mediated IL-10 synthesis",
                location = 0.55)
p4
p_all <- p1 + p2 + p3 + p4 +
  plot_layout(
    ncol = 4
  )

p_all
ggsave(filename = "Fig7C-GSEA.png",
       plot = p_all, width = 16, height = 3, dpi = 300,bg = "white")
ggsave(filename = "Fig7C-GSEA.pdf",
       plot = p_all, width = 16, height = 3, bg = "white",device = cairo_pdf)

#--------7D--------
sender_diff <- read.csv(
  file = "sender_diff.csv"
)

sender_diff$pair <- paste0(
  sender_diff$ligand,
  " -> ",
  sender_diff$receptor
)

sender_pairs <- c(
  "CD80 -> CTLA4",
  "CD80 -> CD274",
  "ICOSLG -> CTLA4",
  "ICOSLG -> ICOS",
  "CXCL9 -> CXCR3",
  "CXCL10 -> CXCR3",
  "CXCL11 -> CXCR3",
  "CCL8 -> CCR1",
  "CCL5 -> CCR5",
  "TNF -> TNFRSF1B",
  "TNF -> TNFRSF1A",
  "IL10 -> IL10RA",
  "IL10 -> IL10RB",
  "IL1B -> IL1R2",
  "GAS6 -> AXL",
  "GAS6 -> MERTK",
  "C3 -> C3AR1",
  "C3 -> ITGAX_ITGB2",
  "THBS1 -> CD47"
)

sender_diff <- sender_diff %>%
  filter(pair %in% sender_pairs)

sender_diff_high <- sender_diff %>%
  select(source,target,pair,prob_high,pathway_name) %>%
  mutate(`log2(prob)` = log2(prob_high)) %>%
  select(-prob_high)

sender_diff_high$source <- "TCN2-high TAM_LYVE1"

sender_diff_low <- sender_diff %>%
  select(source,target,pair,prob_low,pathway_name) %>%
  mutate(`log2(prob)` = log2(prob_low)) %>%
  select(-prob_low)

sender_diff_low$source <- "TCN2-low TAM_LYVE1"

sender_diff <- rbind(
  sender_diff_high,
  sender_diff_low
)

sender_diff$`log2(prob)`[
  sender_diff$`log2(prob)` == -Inf
] <- NA

sender_diff <- na.omit(sender_diff)

celltype_order <- c(
  "Th17",
  "Treg",
  "CD8+ Tex",
  "Tfh",
  "Proliferating T",
  "MAIT",
  "CD8+ Teff",
  "CD4+ Tmem",
  "CD4+ Tnaive",
  "CD8+ Temra",
  "NK_CD16high",
  "NK_CD16low",
  "TCN2-high",
  "TCN2-low",
  "TAM_IL1A",
  "Mono_classical",
  "Mono_nonclassical",
  "Neutrophil",
  "cDC1",
  "cDC2",
  "tDC",
  "pDC",
  "Fibroblast",
  "Smooth Muscle cell",
  "Endothelial",
  "Mast cell",
  "Plasma cell",
  "Epithelial"
)

plot_df <- sender_diff %>%
  mutate(
    col_id = paste(source,target,sep="_")
  )

plot_df$pair <- factor(
  plot_df$pair,
  levels = sender_pairs
)

plot_df$target <- factor(
  plot_df$target,
  levels = celltype_order
)

plot_df$source <- factor(
  plot_df$source,
  levels = c(
    "TCN2-high TAM_LYVE1",
    "TCN2-low TAM_LYVE1"
  )
)

plot_df <- plot_df %>%
  arrange(
    source,
    target
  ) %>%
  mutate(
    col_id = factor(
      col_id,
      levels = unique(col_id)
    )
  )

mat_df <- plot_df %>%
  select(
    pair,
    col_id,
    `log2(prob)`
  ) %>%
  pivot_wider(
    names_from = col_id,
    values_from = `log2(prob)`
  ) %>%
  arrange(
    match(pair,sender_pairs)
  )
mat_df <- as.data.frame(mat_df)
rownames(mat_df) <- mat_df$pair
mat <- as.matrix(mat_df[,-1])

anno_df <- plot_df %>%
  distinct(
    col_id,
    source,
    target
  ) %>%
  arrange(
    match(col_id,colnames(mat))
  )

celltype_colors <- c(
  "Th17"="#6A5ACD",
  "Treg"="#836FFF",
  "CD8+ Tex"="#4169E1",
  "Tfh"="#6495ED",
  "Proliferating T"="#87CEFA",
  "MAIT"="#4682B4",
  "CD8+ Teff"="#1E90FF",
  "CD4+ Tmem"="#5F9EA0",
  "CD4+ Tnaive"="#87CEEB",
  "CD8+ Temra"="#27408B",
  "NK_CD16high"="#20B2AA",
  "NK_CD16low"="#66CDAA",
  "TCN2-high"="#D73027",
  "TCN2-low"="#FC8D59",
  "TAM_IL1A"="#E6550D",
  "Mono_classical"="#F16913",
  "Mono_nonclassical"="#FD8D3C",
  "Neutrophil"="#FDBB84",
  "cDC1"="#B8860B",
  "cDC2"="#DAA520",
  "tDC"="#E6AB02",
  "pDC"="#FDD835",
  "Fibroblast"="#238B45",
  "Smooth Muscle cell"="#41AB5D",
  "Endothelial"="#74C476",
  "Mast cell"="#E78AC3",
  "Plasma cell"="#A6761D",
  "Epithelial"="#7F7F7F"
)

ha <- HeatmapAnnotation(
  Source = anno_df$source,
  Target = anno_df$target,
  col = list(
    Source = c(
      "TCN2-high TAM_LYVE1"="#8B5A9E",
      "TCN2-low TAM_LYVE1"="#A6CEE3"
    ),
    Target = celltype_colors
  ),
  annotation_name_side = "right",
  annotation_name_gp = gpar(
    fontsize=10,
    fontface="bold"
  ),
  annotation_height = unit(
    c(6,6),
    "mm"
  )
)

row_anno_df <- plot_df %>%
  distinct(
    pair,
    pathway_name
  ) %>%
  arrange(
    match(pair,sender_pairs)
  )
pathway_colors <- c(
  "CCL"       = "#8FA8BA",
  "CD80"      = "#365C2D",
  "COMPLEMENT"= "#4B2A63",
  "CXCL"      = "#B2A2BC",
  "GAS"       = "#A8D0E6",
  "ICOS"      = "#C5B6C5",
  "IL1"       = "#806D7A",
  "THBS"      = "#17351E",
  "TNF"       = "#6A2446"
)
left_ha <- rowAnnotation(
  Pathway = row_anno_df$pathway_name,
  col = list(
    Pathway = pathway_colors
  ),
  show_annotation_name = FALSE,
  annotation_width = unit(
    25,
    "mm"
  )
)

heatmap_col <- colorRamp2(
  seq(
    min(mat,na.rm=TRUE),
    max(mat,na.rm=TRUE),
    length.out=6
  ),
  c(
    "#4B65B0",
    "#7FCBAE",
    "#D9E88F",
    "#FDD97D",
    "#F46F43",
    "#A4054E"
  )
)

p <- Heatmap(
  mat,
  name="log2(prob)",
  col=heatmap_col,
  na_col="white",
  top_annotation=ha,
  left_annotation=left_ha,
  cluster_rows=FALSE,
  cluster_columns=FALSE,
  show_column_names=FALSE,
  show_row_names=TRUE,
  row_names_side="right",
  row_names_gp=gpar(
    fontsize=12
  ),
  heatmap_legend_param=list(
    title="log2(prob)",
    title_position="topcenter"
  ),
  border=TRUE
)

pdf(file = "Fig7D-sender.pdf", width = 10, height = 5)
draw(
  p
)
dev.off()

#--------7E--------
receiver_diff <- read.csv(
  file = "receiver_diff.csv"
)

receiver_diff$pair <- paste0(
  receiver_diff$ligand,
  " -> ",
  receiver_diff$receptor
)

receiver_pairs <- c(

  "FN1 -> ITGA4_ITGB7",
  "FN1 -> ITGA5_ITGB1",
  "COL1A1 -> ITGA2_ITGB1",
  "POSTN -> ITGAV_ITGB5",
  "THBS1 -> CD47",

  "GAS6 -> AXL",
  "GAS6 -> MERTK",

  "IL1B -> IL1R2",
  "IL1B -> IL1R1",
  "TNF -> TNFRSF1A",
  "TNF -> TNFRSF1B",

  "CCL3 -> CCR5",
  "CCL4 -> CCR5",
  "CXCL12 -> CXCR4",

  "MIF -> CD74_CXCR4"
)

receiver_diff <- receiver_diff %>%
  filter(pair %in% receiver_pairs)

receiver_diff_high <- receiver_diff %>%
  select(source,target,pair,prob_high,pathway_name) %>%
  mutate(`log2(prob)` = log2(prob_high)) %>%
  select(-prob_high)

receiver_diff_high$target <- "TCN2-high TAM_LYVE1"

receiver_diff_low <- receiver_diff %>%
  select(source,target,pair,prob_low,pathway_name) %>%
  mutate(`log2(prob)` = log2(prob_low)) %>%
  select(-prob_low)

receiver_diff_low$target <- "TCN2-low TAM_LYVE1"

receiver_diff <- rbind(
  receiver_diff_high,
  receiver_diff_low
)

receiver_diff$`log2(prob)`[
  receiver_diff$`log2(prob)` == -Inf
] <- NA

receiver_diff <- na.omit(receiver_diff)

celltype_order <- c(
  "Th17",
  "Treg",
  "CD8+ Tex",
  "Tfh",
  "Proliferating T",
  "MAIT",
  "CD8+ Teff",
  "CD4+ Tmem",
  "CD4+ Tnaive",
  "CD8+ Temra",
  "NK_CD16high",
  "NK_CD16low",
  "TCN2-high",
  "TCN2-low",
  "TAM_IL1A",
  "Mono_classical",
  "Mono_nonclassical",
  "Neutrophil",
  "cDC1",
  "cDC2",
  "tDC",
  "pDC",
  "Fibroblast",
  "Smooth Muscle cell",
  "Endothelial",
  "Mast cell",
  "Plasma cell",
  "Epithelial"
)

plot_df <- receiver_diff %>%
  mutate(
    col_id = paste(source,target,sep="_")
  )

plot_df$pair <- factor(
  plot_df$pair,
  levels = receiver_pairs
)

plot_df$source <- factor(
  plot_df$source,
  levels = celltype_order
)

plot_df$target <- factor(
  plot_df$target,
  levels = c(
    "TCN2-high TAM_LYVE1",
    "TCN2-low TAM_LYVE1"
  )
)

plot_df <- plot_df %>%
  arrange(
    target,
    source
  ) %>%
  mutate(
    col_id = factor(
      col_id,
      levels = unique(col_id)
    )
  )

mat_df <- plot_df %>%
  select(
    pair,
    col_id,
    `log2(prob)`
  ) %>%
  pivot_wider(
    names_from = col_id,
    values_from = `log2(prob)`
  ) %>%
  arrange(
    match(pair,receiver_pairs)
  )
mat_df <- as.data.frame(mat_df)
rownames(mat_df) <- mat_df$pair
mat <- as.matrix(mat_df[,-1])

anno_df <- plot_df %>%
  distinct(
    col_id,
    source,
    target
  ) %>%
  arrange(
    match(col_id,colnames(mat))
  )

celltype_colors <- c(
  "Th17"="#6A5ACD",
  "Treg"="#836FFF",
  "CD8+ Tex"="#4169E1",
  "Tfh"="#6495ED",
  "Proliferating T"="#87CEFA",
  "MAIT"="#4682B4",
  "CD8+ Teff"="#1E90FF",
  "CD4+ Tmem"="#5F9EA0",
  "CD4+ Tnaive"="#87CEEB",
  "CD8+ Temra"="#27408B",
  "NK_CD16high"="#20B2AA",
  "NK_CD16low"="#66CDAA",
  "TCN2-high"="#D73027",
  "TCN2-low"="#FC8D59",
  "TAM_IL1A"="#E6550D",
  "Mono_classical"="#F16913",
  "Mono_nonclassical"="#FD8D3C",
  "Neutrophil"="#FDBB84",
  "cDC1"="#B8860B",
  "cDC2"="#DAA520",
  "tDC"="#E6AB02",
  "pDC"="#FDD835",
  "Fibroblast"="#238B45",
  "Smooth Muscle cell"="#41AB5D",
  "Endothelial"="#74C476",
  "Mast cell"="#E78AC3",
  "Plasma cell"="#A6761D",
  "Epithelial"="#7F7F7F"
)

ha <- HeatmapAnnotation(
  Source = anno_df$source,
  Target = anno_df$target,
  col = list(
    Target = c(
      "TCN2-high TAM_LYVE1"="#8B5A9E",
      "TCN2-low TAM_LYVE1"="#A6CEE3"
    ),
    Source = celltype_colors
  ),
  annotation_name_side = "right",
  annotation_name_gp = gpar(
    fontsize=10,
    fontface="bold"
  ),
  annotation_height = unit(
    c(6,6),
    "mm"
  )
)

row_anno_df <- plot_df %>%
  distinct(
    pair,
    pathway_name
  ) %>%
  arrange(
    match(pair,receiver_pairs)
  )
pathway_colors <- c(
  "CCL"       = "#8FA8BA",
  "CD80"      = "#365C2D",
  "COMPLEMENT"= "#4B2A63",
  "CXCL"      = "#B2A2BC",
  "GAS"       = "#A8D0E6",
  "ICOS"      = "#C5B6C5",
  "IL1"       = "#806D7A",
  "THBS"      = "#17351E",
  "TNF"       = "#6A2446",
  "FN1"        = "#4E79A7",
  "PERIOSTIN"  = "#59A14F",
  "SPP1"       = "#E15759",
  "MIF"        = "#F28E2B"
)
left_ha <- rowAnnotation(
  Pathway = row_anno_df$pathway_name,
  col = list(
    Pathway = pathway_colors
  ),
  show_annotation_name = FALSE,
  annotation_width = unit(
    25,
    "mm"
  )
)

heatmap_col <- colorRamp2(
  seq(
    min(mat,na.rm=TRUE),
    max(mat,na.rm=TRUE),
    length.out=6
  ),
  c(
    "#4B65B0",
    "#7FCBAE",
    "#D9E88F",
    "#FDD97D",
    "#F46F43",
    "#A4054E"
  )
)

p <- Heatmap(
  mat,
  name="log2(prob)",
  col=heatmap_col,
  na_col="white",
  top_annotation=ha,
  left_annotation=left_ha,
  cluster_rows=FALSE,
  cluster_columns=FALSE,
  show_column_names=FALSE,
  show_row_names=TRUE,
  row_names_side="right",
  row_names_gp=gpar(
    fontsize=12
  ),
  heatmap_legend_param=list(
    title="log2(prob)",
    title_position="topcenter"
  ),
  border=TRUE
)

p
pdf(file = "Fig7E-receiver.pdf", width = 10, height = 5)
draw(
  p
)
dev.off()

