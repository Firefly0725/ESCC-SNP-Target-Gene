library(xlsx)
library(showtext)
library(dplyr)
library(circlize)
library(tidyr)
library(ComplexHeatmap)
library(grid)
library(gridExtra)
library(ggplot2)
library(cowplot)

rm(list = ls()); gc()
ORIGINAL_DIR <- ""
setwd(ORIGINAL_DIR)

#--------2A--------
result <- read.table("C:/Users/GengR/Downloads/type1.txt")
result <- result %>%  mutate(
  length = case_when(
    V3 >= V7 ~ V3 - V7,
    V3 <= V6 ~ V6 - V3,
    TRUE ~ 0
  )
)
result <- result %>% select(V1,V4,V9,length)
result <- result %>%
  group_by(V4,V9) %>%
  slice_max(length, n = 1, with_ties = FALSE)


#type1_chr2
pdf("Fig2A_chord_diagram_chr2.pdf", width = 17, height = 15) 
circos.par(gap.after = 5)  
circos.clear()

data1 <- result %>% filter(V1 == "chr2") %>% select(V4, V9, length)
min_len <- min(data1$length, na.rm = TRUE)
max_len <- max(data1$length, na.rm = TRUE)
if (min_len == max_len) {
  data1$length <- 1 
} else {
  data1$length <- 1 + (data1$length - min_len) * 99 / (max_len - min_len)
}
data1$length <- round(data1$length)
colnames(data1) <- c("X4","X8","X11")
data1[c(14, 2), ] <- data1[c(2, 14), ]
data1[c(15, 5), ] <- data1[c(5, 15), ]

grid_col <- c("rs35963873" = "plum1", 
              "rs142741123" = "darkorchid4", 
              "rs3769025" = "blue4", 
              "rs3769823" = "seagreen1")
chordDiagram(data1, annotationTrack = c("grid"),
             grid.col=grid_col,
             order = c("rs3769823", "rs3769025", "rs142741123", "rs35963873", data1$X8))

circos.track(track.index = 1, panel.fun = function(x, y) {
  sector.index <- get.cell.meta.data("sector.index")
  xlim <- get.cell.meta.data("xlim")
  ylim <- get.cell.meta.data("ylim")
  
  if (sector.index=="rs142741123" || sector.index=="rs35963873" || sector.index=="rs3769025" || sector.index=="rs3769823") {
    circos.text(mean(xlim), ylim[2] + mm_y(5), sector.index,
                facing = "outside",  
                niceFacing = TRUE,       
                adj = c(0, 0.5),         
                cex = 1.8,
                col = "black")          
  } else {
    circos.text(mean(xlim), ylim[2] + mm_y(1), sector.index,
                facing = "clockwise",  
                niceFacing = TRUE,            
                adj = c(0, 0.5),            
                cex = 1.0)                
  }
})
circos.clear()
dev.off()

#type1_chr22
pdf("Fig2A_chord_diagram_chr22.pdf", width = 17, height = 15) 
data1 <- result %>% filter(V1 == "chr22") %>% select(V4, V9, length)
min_len <- min(data1$length, na.rm = TRUE)
max_len <- max(data1$length, na.rm = TRUE)
if (min_len == max_len) {
  data1$length <- 1 
} else {
  data1$length <- 1 + (data1$length - min_len) * 99 / (max_len - min_len)
}
data1$length <- round(data1$length)
colnames(data1) <- c("X4","X8","X11")
data1[c(23,3), ] <- data1[c(3, 23), ]
data1[c(25, 62), ] <- data1[c(62, 25), ]
data1[c(26, 16), ] <- data1[c(16, 26), ]

circos.par(gap.after = 5)  
circos.clear()

grid_col <- c("rs10854810" = "royalblue3", 
              "rs3788409" = "aquamarine1", 
              "rs5753220" = "darkgoldenrod3", 
              "rs6005863" = "moccasin")
chordDiagram(data1,annotationTrack = c("grid"),
             grid.col = grid_col)

circos.track(track.index = 1, panel.fun = function(x, y) {
  sector.index <- get.cell.meta.data("sector.index")
  xlim <- get.cell.meta.data("xlim")
  ylim <- get.cell.meta.data("ylim")
  
  if (sector.index=="rs10854810" || sector.index=="rs3788409" || sector.index=="rs5753220" || sector.index=="rs6005863") {
  
    circos.text(mean(xlim), ylim[2] + mm_y(5), sector.index,
                facing = "outside",   
                niceFacing = TRUE,      
                adj = c(0, 0.5),       
                cex = 1.8,
                col = "black")              
  } else {
    circos.text(mean(xlim), ylim[2] + mm_y(1), sector.index,
                facing = "clockwise",  
                niceFacing = TRUE,            
                adj = c(0, 0.5),            
                cex = 1.0)               
  }
})
circos.clear()
dev.off()


#--------2B--------
#data preparation
#type1
cellline <- c("KYSE70","KYSE140","TT")
for(k in cellline){
  data <- read.table(paste0("type1_",k,".txt"))
  data$snpTarget <- paste(data$V4, data$V5, sep = "_"); data$cellLine <- k; data$value <- 1
  data <- data %>% select(snpTarget, cellLine, value)
  assign(paste0("type1_",k),data)
}
df_list <- list(
  KYSE70 = type1_KYSE70,
  KYSE140 = type1_KYSE140,
  TT = type1_TT
)

combined <- bind_rows(df_list, .id = "cellLine") 
result <- combined %>%
  select(snpTarget, cellLine) %>%   
  mutate(exists = 1) %>%                 
  pivot_wider(names_from = cellLine, 
              values_from = exists, 
              values_fill = 0) 
result$KYSE180 <- 0
result$type <- "type1"
result$number <- result$KYSE140 + result$KYSE180 + result$KYSE70 + result$TT
result_type1 <- result %>% select(KYSE140,KYSE180, KYSE70, TT)
anno_type1 <- result %>% select(type, number)

#type2
cellline <- c("KYSE70","KYSE140","TT")
for(k in cellline){
  data <- read.table(paste0("type2_",k,".txt"))
  data$snpTarget <- paste(data$V4, data$V5, sep = "_"); data$cellLine <- k; data$value <- 1
  data <- data %>% select(snpTarget, cellLine, value)
  assign(paste0("type2_",k),data)
}
df_list <- list(
  KYSE70 = type2_KYSE70,
  KYSE140 = type2_KYSE140,
  TT = type2_TT
)

combined <- bind_rows(df_list, .id = "cellLine") 
result <- combined %>%
  select(snpTarget, cellLine) %>%   
  mutate(exists = 1) %>%                 
  pivot_wider(names_from = cellLine, 
              values_from = exists, 
              values_fill = 0) 
result$KYSE180 <- 0
result$type <- "type2"
result$number <- result$KYSE140 + result$KYSE180 + result$KYSE70 + result$TT
result_type2 <- result %>% select(KYSE140,KYSE180, KYSE70, TT)
anno_type2 <- result %>% select(type, number)

#type3
cellline <- c("KYSE70","KYSE140","TT","KYSE180")

for(k in cellline){
  data <- read.table(paste0("type3_",k,".txt"))
  data$snpTarget <- paste(data$V4, data$V5, sep = "_")
  data$cellLine <- k
  data$value <- 1
  data <- data %>% select(snpTarget, cellLine, value)
  assign(paste0("type3_",k), data)
}
combined <- bind_rows(type3_KYSE70, type3_KYSE140, type3_TT, type3_KYSE180)
combined <- unique(combined)
result <- combined %>%
  select(snpTarget, cellLine) %>%   
  mutate(exists = 1) %>%                 
  pivot_wider(names_from = cellLine, 
              values_from = exists, 
              values_fill = 0) 
result$type <- "type3"
result$number <- result$KYSE140 + result$KYSE180 + result$KYSE70 + result$TT
result_type3 <- result %>% select(KYSE140,KYSE180, KYSE70, TT)
anno_type3 <- result %>% select(type, number)

#figure
data <- rbind(result_type1,result_type2)
col_sums <- rowSums(result_type3, na.rm = TRUE) 
counts <- table(factor(col_sums, levels = 1:4))


data <- t(data)
data <- as.matrix(data)
anno <- rbind(anno_type1, anno_type2)
anno <- t(anno)

top_anno = HeatmapAnnotation(type = anno[1,], number = anno[2,], border = T,
                             show_annotation_name = F,
                             col = list(
                               type = c("type1" = "darkolivegreen1", "type2" = "mediumturquoise"),
                               number = c("1" = "lightblue1", "2" = "turquoise1", "3" = "bisque", "4" = "darkorchid3")
                             ),
                             annotation_legend_param = list(
                               type = list(title = "Types", 
                                           title_gp = gpar(fontsize = 18, fontface = "bold", col = "black"),
                                           labels_gp = gpar(fontsize = 18, col = "black")),
                               number = list(title = "Numbers", 
                                             title_gp = gpar(fontsize = 18, fontface = "bold", col = "black"),
                                             labels_gp = gpar(fontsize = 18, col = "black"))
                               
                             ))
p1 <- Heatmap(data,cluster_rows = F,
              cluster_columns = F,
              show_column_names = F,
              show_column_dend = F,
              col = c('white','deeppink2'),
              border = 'black',
              show_heatmap_legend = F,
              top_annotation =top_anno,
              row_names_side = "left",  
              row_names_gp = gpar(fontsize = 18))


data1 <- result_type3
data1 <- t(data1)
data1 <- as.matrix(data1)
anno1 <- anno_type3
anno1 <- t(anno1)

top_anno = HeatmapAnnotation(type = anno1[1,], number = anno1[2,], border = T,
                             show_annotation_name = F,
                             col = list(
                               type = c("type3" = "deepskyblue"),
                               number = c("1" = "lightblue1", "2" = "turquoise1", "3" = "bisque", "4" = "darkorchid3")
                             ),
                             annotation_legend_param = list(
                               type = list(title = "Types", 
                                           title_gp = gpar(fontsize = 18, fontface = "bold", col = "black"),
                                           labels_gp = gpar(fontsize = 18, col = "black")),
                               number = list(title = "Numbers", 
                                             title_gp = gpar(fontsize = 18, fontface = "bold", col = "black"),
                                             labels_gp = gpar(fontsize = 18, col = "black"))
                               
                             ))
p2 <- Heatmap(data1,cluster_rows = F,
              cluster_columns = F,
              show_column_names = F,
              show_column_dend = F,
              col = c('white','deeppink2'),
              border = 'black',
              show_heatmap_legend = F,
              top_annotation =top_anno,
              row_names_side = "left",    
              row_names_gp = gpar(fontsize = 18))
p2 

g1 <- grid.grabExpr(draw(p1))
g2 <- grid.grabExpr(draw(p2))
p <- grid.arrange(g1, g2, ncol = 1)
ggsave(filename = file.path("Fig2B-snp_target.png"), 
       plot = p, width = 10, height = 5, dpi = 300,bg = "white")
ggsave(filename = file.path("Fig2B-snp_target.pdf"), 
       plot = p, width = 10, height = 5,bg = "white")


#--------2C--------
data <- data.frame(
  Category = c("type 1", "type 1", "type 3", "type 3"),
  Group    = c("cis", "trans", "cis", "trans"),
  Value    = c(10, 2, 17, 6)
)

data <- data %>%
  group_by(Category) %>%
  mutate(
    total = sum(Value),
    prop = Value / total,
    label = sprintf("%.1f%% (%d)", prop * 100, Value)
  ) %>%
  ungroup()

plot1 <- ggplot(data, aes(x = Category, y = Value, fill = Group)) +
  geom_bar(stat = "identity", position = "fill", colour = "black") +
  geom_text(aes(label = label), 
            position = position_fill(vjust = 0.5), 
            size = 14 / .pt,      
            colour = "black") +
  xlab("Types") + 
  ylab("Proportion") +
  scale_fill_manual(values = c("cis" = "#87CEEB", 
                               "trans" = "#E64B35CC"),
                    name = "eQTL type") +
  scale_y_continuous(labels = scales::percent) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        panel.border = element_rect(color = "black", fill = NA, size = 1),
        axis.text.x = element_text(size = 14,color = "black"),
        axis.text.y = element_text(size = 14,color = "black"),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 14),
        axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14))
plot1
ggsave(filename = file.path("Fig2C-eQTL_result.png"), 
       plot = plot1, width = 5, height = 5, dpi = 300,bg = "white")
ggsave(filename = file.path("Fig2C-eQTL_result.pdf"), 
       plot = plot1, width = 5, height = 5, bg = "white")
