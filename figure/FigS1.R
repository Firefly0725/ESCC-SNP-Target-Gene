rm(list = ls()); gc()
ORIGINAL_DIR <- ""
setwd(ORIGINAL_DIR)
output <- file.path("")
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


#----S1A----
Myeloid_sc_data <- readRDS("merge_seurat_Myeloid.RDS")
p1 <- DimPlot(Myeloid_sc_data, group.by = "cell_type")+
  ggtitle("UMAP of Celltype") +             
  theme(
    plot.title = element_text(size = 16),   
    axis.title = element_text(size = 16), 
    axis.text = element_text(size = 16),
    legend.text = element_text(size = 15)
  )
ggsave(file.path("Myeloid_UMAP.png"), p1, width = 7, height = 5, dpi = 300, bg = "white")
ggsave(file.path("Myeloid_UMAP.pdf"), p1, width = 7, height = 5, bg = "white")

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
    axis.text.x = element_text(size = 14, angle = 45),  
    axis.text.y = element_text(size = 16),             
    axis.title = element_text(size = 16),
    legend.text = element_text(size = 14),
    legend.title = element_text(size = 14),
    strip.text.x = element_text(size = 14)
  ) +
  xlab(NULL) + 
  ylab("Cell Type")
ggsave(file.path("Myeloid_dotplot.png"), p1, width = 11, height = 4, dpi = 300, bg = "white")
ggsave(file.path("Myeloid_dotplot.pdf"), p1, width = 11, height = 4, bg = "white")


#--------S1B--------
CD4_sc_data <- readRDS("merge_seurat_Tcell.RDS")
p1 <- DimPlot(CD4_sc_data, group.by = "cell_type",label.size = 6)+
  ggtitle("UMAP of celltype") +             
  theme(
    plot.title = element_text(size = 16),   
    axis.title = element_text(size = 16), 
    axis.text = element_text(size = 16),
    legend.text = element_text(size = 15)
  )
ggsave(file.path("Tcell_UMAP.png"), p1, width = 7, height = 5, dpi = 300, bg = "white")
ggsave(file.path("Tcell_UMAP.pdf"), p1, width = 7, height = 5, bg = "white")

marker_genes <- c(
  "CCR7","SELL","LEF1","S1PR1",
  "LMNA","ANXA1","FOS",
  "CD4","CCR6","KLRB1","CTSH",
  "CTSL","PTPN13","TOX2","NMB","PASK","GNG4",
  "FOXP3","CTLA4","IL2RA",
  "CD8A","CD8B","ISG20",
  "CX3CR1","HOPX","ZEB2",
  "CXCL13","TOX","HAVCR2","LAYN","CCL3","VCAM1","RGS2",
  "SLC4A10","ZBTB16","RORA","RORC",
  "MKI67","TOP2A"
)
p1 <- DotPlot(CD4_sc_data, features = marker_genes, group.by = "cell_type")+
  RotatedAxis()+    
  scale_color_gradientn(colors = c("blue", "white", "red"))+
  theme(
    axis.text.x = element_text(size = 14, angle = 45),
    axis.text.y = element_text(size = 16),            
    axis.title = element_text(size = 16),
    legend.text = element_text(size = 14),
    legend.title = element_text(size = 14),
    strip.text.x = element_text(size = 14)
  ) +
  xlab(NULL) + 
  ylab("Cell Type")
ggsave(file.path("Tcell_dotplot.png"), p1, width = 12, height = 4, dpi = 300, bg = "white")
ggsave(file.path("Tcell_dotplot.pdf"), p1, width = 12, height = 4, bg = "white")


#--------S1C--------
CD4_sc_data <- readRDS("merge_seurat_NK.RDS")
p1 <- DimPlot(CD4_sc_data, group.by = "cell_type",label.size = 6)+
  ggtitle("UMAP of celltype") +             
  theme(
    plot.title = element_text(size = 16),   
    axis.title = element_text(size = 16), 
    axis.text = element_text(size = 16),
    legend.text = element_text(size = 15)
  )
ggsave(file.path("NK_UMAP.png"), p1, width = 6, height = 5, dpi = 300, bg = "white")
ggsave(file.path("NK_UMAP.pdf"), p1, width = 6, height = 5, bg = "white")

marker_genes <- c("FCGR3A","PRF1","GZMB",
                  "NCAM1","GZMK", "TCF7", "SELL")
p1 <- DotPlot(CD4_sc_data, features = marker_genes, group.by = "cell_type")+
  RotatedAxis()+    
  scale_color_gradientn(colors = c("blue", "white", "red"))+
  theme(
    axis.text.x = element_text(size = 14, angle = 45), 
    axis.text.y = element_text(size = 16),            
    axis.title = element_text(size = 16),
    legend.text = element_text(size = 14),
    legend.title = element_text(size = 14),
    strip.text.x = element_text(size = 14)
  ) +
  xlab(NULL) + 
  ylab("Cell Type")
ggsave(file.path("NK_dotplot.png"), p1, width = 7, height = 4, dpi = 300, bg = "white")
ggsave(file.path("NK_dotplot.pdf"), p1, width = 7, height = 4, bg = "white")


