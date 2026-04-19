rm(list = ls()); gc()
ORIGINAL_DIR <- "D:/大学的资料/R/ESCC/data/scRNA/"
setwd(ORIGINAL_DIR)
output <- file.path("D:/大学的资料/R/ESCC/data/scRNA/")
if (!dir.exists(output)) dir.create(output, recursive = TRUE)
library(Seurat)
library(ggplot2)
library(dplyr)
library(harmony)
library(stringr)
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(tidyr)
library(pheatmap) 
library(ggpubr)
filtered_seurat_umap <- readRDS(file = paste0(output,"/filtered_seurat_umap.rds"))

target_cells <- subset(filtered_seurat_umap, subset = cell_type == "Myeloid")
matrix <- GetAssayData(target_cells, assay = "RNA", layer = "count")
Myeloid_sc_data <- CreateSeuratObject(
  counts = matrix,  
  project = "scRNA_project", 
  min.cells = 3,            
  min.features = 200,
  meta.data = target_cells@meta.data
)
Myeloid_sc_data@meta.data <- Myeloid_sc_data@meta.data[, !(colnames(Myeloid_sc_data@meta.data) %in% 
                                                             c("cluster", "cell_type","seurat_clusters","RNA_snn_res.0.6"))]


Myeloid_sc_data <- NormalizeData(Myeloid_sc_data, normalization.method = "LogNormalize", scale.factor = 1e4)
Myeloid_sc_data <- FindVariableFeatures(Myeloid_sc_data, selection.method = "vst")
Myeloid_sc_data[["RNA"]] <- subset(Myeloid_sc_data[["RNA"]], layers = "scale.data")
Myeloid_sc_data <- ScaleData(Myeloid_sc_data)
Myeloid_sc_data <- RunPCA(Myeloid_sc_data,npcs = 50)
ElbowPlot(Myeloid_sc_data, ndims = 50)

Myeloid_sc_data <- RunHarmony(Myeloid_sc_data, group.by.vars = "Sample")
Myeloid_sc_data <- FindNeighbors(Myeloid_sc_data,reduction = "harmony",dims = 1:35)
Myeloid_sc_data <- FindClusters(Myeloid_sc_data, resolution = 1.2)
Myeloid_sc_data <- RunUMAP(Myeloid_sc_data, dims = 1:35, reduction = "harmony",
                           n.neighbors = 10, min.dist = 0.5, spread =0.6)
Myeloid_sc_data$cluster <- Idents(Myeloid_sc_data) 
DimPlot(Myeloid_sc_data, group.by = "cluster",label = TRUE, label.size = 6)+
  ggtitle("UMAP of Samples") +             
  theme(
    plot.title = element_text(size = 16),   
    axis.title = element_text(size = 16), 
    axis.text = element_text(size = 16),
    legend.text = element_text(size = 15)
  )
marker_genes <- list(
  "Mono_classical"=c("FCN1","CD14","CD36","SELL","CLEC12A"),
  "Mono_intermediate"=c("FCGR3A","LST1","LILRB2","CD300E","FAM65B"),
  "Mono_nonclassical"=c("FCGR3B","LILRA","LILRB1","CX3XR1","IFITM1","ICAM2"),
  "Macro_RTM"=c("CD80","CD86","NLRP3","LYVE1"),
  #"Macro_M2"=c("CD206","EREG","SPP1","C1QC","CD163","FPR3","C1QA","CCL18","IL10"),
  "Macro_TAM"=c("C1QC","MERTK","FPR3","TREM2","SPP1","IL1RN","OLR1","EREG","VEGFA"),
  "M2-like"=c("CD163","MRC1","CCL18","IL10"),
  "cDC1"=c("CLEC9A","BATF3","XCR1","IRF8"),
  "cDC2"=c("CD1C","FCER1A","FCER2B","CLEC10A"),
  "tDC"=c("FLT3","IDO1","PD-L1","PD-L2"),
  "pDC"=c("LILRA4","GZMB","IL3RA","SLC15A4"))

marker_genes <- list(
  "M1"=c("CCL5","CCR7","CD40","CD86","CXCL9","CXCL10","CXCL11","IDO1","IL1A","IL1B","IL6","IRF1","IRF5","KYNU"),
  "M2"=c("CCL4","CCL13","CCL18","CCL20","CCL22","CD276","CLEC7A","CTSA","CTSB","CTSC","CTSD",
         "FN1","IL4R","IRF4","LYVE1","MMP9","MMP14","MMP19","MSR1","TGFB1","TGFB2","TGFB3","TNFSF8",
         "TNFSF12","VEGFA","VEGFB","VEGFC","MERTK")
)

DotPlot(Myeloid_sc_data, features = marker_genes, group.by = "cluster")+
  RotatedAxis()+    
  scale_color_gradientn(colors = c("blue", "white", "red"))+
  theme(
    axis.text.x = element_text(size = 14, angle = 45),   # 横轴字体
    axis.text.y = element_text(size = 16),              # 纵轴字体
    axis.title = element_text(size = 16),
    legend.text = element_text(size = 14),
    legend.title = element_text(size = 14),
    strip.text.x = element_text(size = 14)
  ) +
  xlab("Marker Genes") + 
  ylab("Cluster")
FeaturePlot(Myeloid_sc_data, features = "TCN2", pt.size = 0.6, order = TRUE) +
  labs(color = "TCN2") +  
  theme(
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 14),
    legend.text = element_text(size = 14),
    plot.title = element_blank()
  )
DimPlot(Myeloid_sc_data, group.by = "condition",label = TRUE, label.size = 6)+
  ggtitle("UMAP of Samples") +             
  theme(
    plot.title = element_text(size = 16),   
    axis.title = element_text(size = 16), 
    axis.text = element_text(size = 16),
    legend.text = element_text(size = 15)
  )

Myeloid_sc_data <- subset(Myeloid_sc_data, subset = seurat_clusters != 12)

table(Idents(Myeloid_sc_data))
new.cluster.ids <- c("Mono_classical","Mono_classical","TAM_LYVE1","TAM_LYVE1","Mono_classical","TAM_CHIT1",
                     "TAM_LYVE1","TAM_IL1A","Mono_nonclassical","pDC","cDC2",
                     "TAM_IL1A","Mono_classical","tDC","cDC1",
                     "Mono_classical","TAM_IL1A","TAM_IL1A","pDC","cDC2",
                     "tDC","TAM_LYVE1","Mono_classical")
names(new.cluster.ids) <- levels(Myeloid_sc_data)
Myeloid_sc_data <- RenameIdents(Myeloid_sc_data, new.cluster.ids)
Myeloid_sc_data$cell_type <- Idents(Myeloid_sc_data)
Myeloid_sc_data$cell_type <- factor(Myeloid_sc_data$cell_type, 
                                         levels = c("Mono_classical","Mono_intermediate","Mono_nonclassical",
                                                    "TAM_IL1A","TAM_LYVE1","TAM_CHIT1",
                                                    "cDC1","cDC2","pDC","tDC"))
saveRDS(Myeloid_sc_data, "merge_seurat_Myeloid.RDS")

Myeloid_sc_data <- readRDS("merge_seurat_Myeloid.RDS")
#绘制图象
#----p1----
p1 <- DimPlot(Myeloid_sc_data, group.by = "cell_type")+
  ggtitle("UMAP of Celltype") +             
  theme(
    plot.title = element_text(size = 16),   
    axis.title = element_text(size = 16), 
    axis.text = element_text(size = 16),
    legend.text = element_text(size = 15)
  )
ggsave(file.path("D:/大学的资料/R/ESCC/data/FigS1/Myeloid_UMAP.png"), p1, width = 7, height = 5, dpi = 300, bg = "white")
ggsave(file.path("D:/大学的资料/R/ESCC/data/FigS1/Myeloid_UMAP.pdf"), p1, width = 7, height = 5, bg = "white")

#----p2----
marker_genes <- c(
  "FCN1","CD14","CD36","CLEC12A",
  "FCGR3A","LST1","LILRB2","IFITM1","ICAM2",
  "IL1RN","OLR1","EREG","VEGFA","IL1A","FPR3","MERTK","LYVE1","TREM2","C1QC","CHIT1",
  "CLEC9A","BATF3","XCR1",
  "CD1C","FCER1A","CLEC10A",
  "LILRA4","GZMB","IL3RA","SLC15A4",
  "FLT3","IDO1"
  )

p1 <- DotPlot(Myeloid_sc_data, features = marker_genes, group.by = "cell_type")+
  RotatedAxis()+    
  scale_color_gradientn(colors = c("blue", "white", "red"))+
  theme(
    axis.text.x = element_text(size = 14, angle = 45),   # 横轴字体
    axis.text.y = element_text(size = 16),              # 纵轴字体
    axis.title = element_text(size = 16),
    legend.text = element_text(size = 14),
    legend.title = element_text(size = 14),
    strip.text.x = element_text(size = 14)
  ) +
  xlab(NULL) + 
  ylab("Cell Type")
ggsave(file.path("D:/大学的资料/R/ESCC/data/FigS1/Myeloid_dotplot.png"), p1, width = 11, height = 4, dpi = 300, bg = "white")
ggsave(file.path("D:/大学的资料/R/ESCC/data/FigS1/Myeloid_dotplot.pdf"), p1, width = 11, height = 4, bg = "white")

#----p3----
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
ggsave(file.path(output, "Fig7D_TCN2_exp_myeloid.png"), p, width = 6, height = 7, dpi = 300, bg = "white")



p <- FeaturePlot(Myeloid_sc_data, features = "CCL20", pt.size = 0.6, order = TRUE) +
  theme(
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 14),
    legend.text = element_text(size = 14),
    plot.title = element_blank()
  )
p
