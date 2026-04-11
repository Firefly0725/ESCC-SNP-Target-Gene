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

rm(list = ls()); gc()
ORIGINAL_DIR <- ""
setwd(ORIGINAL_DIR)


#--------6A--------
output <- file.path(" ")
filtered_seurat_umap <- readRDS(file = paste0(output,"/filtered_seurat_umap.rds"))
my_colors <- c(
  "Epithelial" = "#0099CC",
  "Fibroblast" = "#99CC99",
  "Smooth Muscle cell" = "#99FFCC",
  "Endothelial" = "#66CCCC",
  "Myeloid" = "#9999FF",
  "T cell" = "#99CCFF",
  "NK cell" = "#CC99CC",
  "B cell" = "#CC3300",
  "Plasma cell" = "#FFCC00",
  "Neutrophil" = "#FF6600",
  "Mast cell" = "#FFCCCC"
)

p1 <- DimPlot(filtered_seurat_umap, 
              group.by = "cell_type", 
              cols = my_colors,
              label = TRUE, 
              label.size = 5) +
  theme(
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 14),
    legend.text = element_text(size = 14),
    plot.title = element_blank(),
    legend.position = "bottom",       
    legend.direction = "horizontal"
  )+
  guides(color = guide_legend(
    ncol = 3,
    override.aes = list(size = 5)     
  ))
p1
ggsave(file.path("Fig6A_celltype_UMAP.png"), p1, width = 6, height = 6.5, dpi = 300, bg = "white")
ggsave(file.path("Fig6A_celltype_UMAP.pdf"), p1, width = 6, height = 6.5, bg = "white")


#--------6B--------
marker_genes <- c(
  "PERP","EPCAM","KRT15","TP63","KRT17",
  "FN1","DCN","COL1A1",
  "MYH11","SORBS1","DMD","KCNMA1",
  "VWF","PLVAP","PECAM1","ENG",
  "CD68","LYZ","CD14","CD163",
  "IL7R","CCR7","FOXP3","CD8A","CD8B","GZMA","GZMH","GNLY",
  "FCGR3A","KLRF1",
  "CD79A","CD79B","MS4A1",
  "CD38","MZB1","IGHG1","JCHAIN","SDC1",
  "FCGR3B","G0S2",
  "KIT","CPA3","TPSAB1"
)
filtered_seurat_umap$Celltype <- filtered_seurat_umap$cell_type
filtered_seurat_umap$Condition <- filtered_seurat_umap$condition
png(file.path("Fig6B-gene_expression.png"), width = 1900, height = 2200, res = 300)
plot_heatmap(dataset = filtered_seurat_umap, 
             markers = marker_genes,
             sort_var = c("Celltype","Sample","Condition"),
             anno_var = c("Celltype","Sample","Condition"),
             anno_colors = list(c("#0099CC","#99CC99","#99FFCC","#66CCCC","#9999FF",
                                  "#99CCFF","#CC99CC","#CC3300","#FFCC00","#FF6600","#FFCCCC"),                                           
                                "Set2",  
                                c("tomato","#009CB8")))
dev.off()

pdf(file.path("Fig6B-gene_expression.pdf"), width = 6, height = 7)
plot_heatmap(dataset = filtered_seurat_umap, 
             markers = marker_genes,
             sort_var = c("Celltype","Sample","Condition"),
             anno_var = c("Celltype","Sample","Condition"),
             anno_colors = list(c("#0099CC","#99CC99","#99FFCC","#66CCCC","#9999FF",
                                  "#99CCFF","#CC99CC","#CC3300","#FFCC00","#FF6600","#FFCCCC"),                                           
                                "Set2",  
                                c("tomato","#009CB8")))
dev.off()


#--------6C--------
p <- FeaturePlot(filtered_seurat_umap, features = "TCN2", pt.size = 0.6, order = TRUE) +
  labs(color = "TCN2") +  
  theme(
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 14),
    legend.text = element_text(size = 14),
    plot.title = element_blank()
  )

ggsave(file.path("Fig6C_TCN2_expression.png"),p, width = 3.5, height = 3, dpi = 300, bg = "white")
ggsave(file.path("Fig6C_TCN2_expression.pdf"),p, width = 3.5, height = 3, bg = "white")


#--------6D--------
umap_coords <- Embeddings(filtered_seurat_umap, reduction = "umap")
df_plot <- as.data.frame(umap_coords)
colnames(df_plot) <- c("UMAP_1", "UMAP_2")
df_plot$condition <- filtered_seurat_umap$condition

df_plot <- df_plot[order(df_plot$condition, decreasing = T), ]

my_colors <- c(
  "Tumor" = "tomato",
  "Normal" = "#009CB8"
)

p1 <- ggplot(df_plot, aes(x = UMAP_1, y = UMAP_2, color = condition)) +
  geom_point(size = 0.03, alpha = 0.8) +
  scale_color_manual(values = my_colors, breaks = c("Tumor", "Normal")) +
  theme_minimal() +
  theme(
    panel.grid = element_blank(),
    axis.line = element_line(color = "black", linewidth = 0.5),
    axis.title = element_text(size = 14, color = "black"),
    axis.text = element_text(size = 14, color = "black"),
    legend.text = element_text(size = 14),
    legend.title = element_text(size = 14),
    plot.title = element_blank(),
    legend.position = "bottom"
    # 删除了 legend.direction = "horizontal"
  ) +
  guides(color = guide_legend(
    override.aes = list(size = 3),
    nrow = 2           # 关键：强制分为两行
  )) +
  labs(x = "umap_1", y = "umap_2")

p1
ggsave(file.path("Fig6D_condition_UMAP.png"), p1, width = 2.7, height = 3.5, dpi = 300, bg = "white")



#--------6E--------
output <- file.path(" ")
Myeloid_sc_data <- readRDS(file = file.path(output,"merge_seurat_Myeloid.RDS"))

expr_data <- FetchData(Myeloid_sc_data, vars = c("TCN2", "cell_type", "condition"))
expr_data <- expr_data[!is.na(expr_data$condition), ]

summary_stats <- expr_data %>%
  group_by(cell_type, condition) %>%
  summarise(
    mean_exp = mean(TCN2),
    sd_exp = sd(TCN2),
    n = n(),
    se_exp = sd_exp / sqrt(n),
    .groups = "drop"
  )
p_values <- expr_data %>%
  group_by(cell_type) %>%
  summarise(
    p_val = wilcox.test(TCN2 ~ condition)$p.value,
    .groups = "drop"
  ) %>%
  mutate(
    p_label = case_when(
      p_val < 0.001 ~ "***",
      p_val < 0.01 ~ "**",
      p_val < 0.05 ~ "*",
      TRUE ~ "ns"
    ),
    y_pos = max(summary_stats$mean_exp + summary_stats$se_exp, na.rm = TRUE) + 0.03
  )

p <- ggplot(summary_stats, aes(x = cell_type, y = mean_exp, fill = condition)) +
  geom_bar(stat = "identity", position = position_dodge(0.9), width = 0.7) +
  geom_errorbar(aes(ymin = mean_exp - se_exp, ymax = mean_exp + se_exp),
                position = position_dodge(0.9), width = 0.2) +
  labs(y = "TCN2 Expression Level", fill = "Condition", x = NULL) +
  theme_minimal() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, size = 1),
    axis.title = element_text(size = 14, color = "black"),
    axis.text.x = element_text(angle = 60, hjust = 1, size = 14, color = "black"),
    axis.text.y = element_text(size = 14, color = "black"),
    legend.text = element_text(size = 14, color = "black"),
    legend.title = element_text(size = 14, color = "black")
  ) +
  scale_fill_manual(values = c("Tumor" = "tomato", "Normal" = "#009CB8"))  

p <- p + geom_text(data = p_values, 
                   aes(x = cell_type, y = y_pos, label = p_label),
                   inherit.aes = FALSE, size = 5, vjust = -0.5)

print(p)
ggsave(file.path("Fig6E_TCN2_exp_myeloid.png"), p, width = 6, height = 7, dpi = 300, bg = "white")
ggsave(file.path("Fig6E_TCN2_exp_myeloid.pdf"), p, width = 6, height = 7,  bg = "white")
                       

#--------6F--------
marker_genes <- list(
  "M1"=c("CCL5","CCR7","CD40","CD86","CXCL9","CXCL10","CXCL11","IDO1","IL1A","IL1B","IL6","IRF1","IRF5","KYNU"),
  "M2"=c("CCL4","CCL13","CCL18","CCL20","CCL22","CD276","CLEC7A","CTSA","CTSB","CTSC","CTSD",
         "FN1","IL4R","IRF4","LYVE1","MMP9","MMP14","MMP19","MSR1","TGFB1","TGFB2","TGFB3","TNFSF8",
         "TNFSF12","VEGFA","VEGFB","VEGFC","MERTK")
)  

tam_cells <- subset(Myeloid_sc_data, subset = cell_type == "TAM_LYVE1")  
genes_to_plot <- marker_genes$M1
genes_to_plot <- genes_to_plot[genes_to_plot %in% rownames(tam_cells)]
expr_matrix <- as.matrix(GetAssayData(tam_cells, slot = "data")[genes_to_plot, , drop = FALSE])
df_plot <- expr_matrix %>%
  as.data.frame() %>%
  tibble::rownames_to_column("gene") %>%
  pivot_longer(cols = -gene, names_to = "cell", values_to = "expression")
p1 <- ggplot(df_plot, aes(x = gene, y = expression, fill = gene)) +
  geom_violin(scale = "width", trim = FALSE, alpha = 0.7) +
  stat_summary(fun = median, geom = "point", size = 1, color = "black") +
  theme_classic() +
  theme(
    axis.title = element_text(size = 14, color = "black"),
    axis.text = element_text(size = 14, color = "black"),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 14, color = "black"), 
    legend.position = "none"
  ) +
  labs(x = "M1 signature", y = "Gene expression")

genes_to_plot <- marker_genes$M2
genes_to_plot <- genes_to_plot[genes_to_plot %in% rownames(tam_cells)]
expr_matrix <- as.matrix(GetAssayData(tam_cells, slot = "data")[genes_to_plot, , drop = FALSE])
df_plot <- expr_matrix %>%
  as.data.frame() %>%
  tibble::rownames_to_column("gene") %>%
  pivot_longer(cols = -gene, names_to = "cell", values_to = "expression")
p2 <- ggplot(df_plot, aes(x = gene, y = expression, fill = gene)) +
  geom_violin(scale = "width", trim = FALSE, alpha = 0.7) +
  stat_summary(fun = median, geom = "point", size = 1, color = "black") +
  theme_classic() +
  theme(
    axis.title = element_text(size = 14, color = "black"),
    axis.text = element_text(size = 14, color = "black"),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 14, color = "black"),
    legend.position = "none"
  ) +
  labs(x = "M2 signature", y = "Gene expression")  
library(patchwork)
combined_plot <- p1 / p2  
print(combined_plot)
ggsave(file.path("Fig6F_M1_M2_sig.png"), combined_plot, width = 9, height = 7, dpi = 300, bg = "white")
ggsave(file.path("Fig6F_M1_M2_sig.pdf"), combined_plot, width = 9, height = 7, bg = "white")

