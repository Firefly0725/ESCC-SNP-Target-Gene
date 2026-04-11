library(ggplot2)
library(patchwork)
library(xlsx)
library(ComplexHeatmap)
library(dplyr)

rm(list = ls()); gc()
ORIGINAL_DIR <- ""
setwd(ORIGINAL_DIR)

#--------3A--------
data <- data.frame(
  category = c("Corresponding to 1 ESCC SNP", "Corresponding to 2 ESCC SNPs", "As an ESCC SNP"),
  value = c(125, 5, 2)
)
data$category <- factor(data$category, 
                        levels = c("Corresponding to 1 ESCC SNP", 
                                   "Corresponding to 2 ESCC SNPs",
                                   "As an ESCC SNP"))

p <- ggplot(data, aes(x = "", y = value, fill = category)) +
  geom_bar(stat = "identity", width = 1) +  # 添加扇形黑色边框
  coord_polar(theta = "y") +
  theme_void() +
  geom_text(aes(label = value), 
            position = position_stack(vjust = 0.5),
            size = 5.5) +  
  scale_fill_manual(values = c("Corresponding to 1 ESCC SNP" = "#F08080", 
                               "Corresponding to 2 ESCC SNPs" = "#ADD8E6", 
                               "As an ESCC SNP" = "#F0E68C"),
                    labels = c("Corresponding to\n1 ESCC SNP",  
                               "Corresponding to\n2 ESCC SNPs",
                               "As an ESCC SNP")) +
  theme(plot.title = element_text(hjust = 0.5), 
        legend.title = element_blank(),
        legend.position = "right",         
        legend.key.height = unit(1.5, "cm"),
        legend.text = element_text(size = 14)
  )
p
ggsave(filename = file.path("Fig3A-pie_1.png"), 
       plot = p, width = 6, height = 3.5, dpi = 300,bg = "white")
ggsave(filename = file.path("Fig3A-pie_1.pdf"), 
       plot = p, width = 6, height = 3.5, bg = "white")


#--------3C--------
data2 <- read.csv(file='heat2.CSV',
                  header=TRUE,sep=',')
data3 <- read.csv(file='heat3.csv',
                  header=TRUE,sep=',')

p1 <- ggplot(data2, aes(cellline,factor))+
  geom_point(aes(size = effectsize,fill=P),shape=21,stroke=0.6) +
  scale_size_continuous(name = 'Effect Size',
                        limit = c(-0.001, 70), 
                        breaks = c(0, 10, 20, 50),
                        range = c(3, 22)) +
  scale_fill_gradient(name = 'P',
                      limit = c(-0.00001,0.050001),
                      breaks = c(0.001,0.01,0.05),
                      low = 'deeppink',
                      high = 'white')+
  scale_y_discrete(limits=as.character(c("SOX2","H3K4me3","H3K4me1","H3K27ac","H3K27ac anchors")))+ 
  scale_x_discrete(limits=as.character(c("KYSE140","KYSE180","KYSE70","TT","TE5")))+ 
  theme_bw()+
  xlab(NULL) + 
  ylab(NULL)+
  theme(panel.border = element_rect(fill=NA,color="black", size=0.5, linetype="solid"),
        axis.text.x = element_text(size = 14, angle = 45, hjust = 1,colour = "black"),  
        axis.text.y = element_text(size = 14, colour = "black"),  
        legend.title = element_text(size = 14),  
        legend.text = element_text(size = 12))  
p1
p2<-p1+geom_point(data=data3,
                  mapping =aes(cellline,factor),
                  shape=4,
                  stroke=1,
                  size=10,
                  color='black')
p2
ggsave(filename = file.path("Fig3C-shuffle.png"), 
       plot = p2, width = 6.5, height = 5, dpi = 300,bg = "white")
ggsave(filename = file.path("Fig3C-shuffle.pdf"), 
       plot = p2, width = 6.5, height = 5, bg = "white")


#--------3D--------
#type1
data <- data.frame(
  category = c("Shared interaction relationships with ESCC SNPs", "Unique interaction relationships"),
  value = c(142, 847)
)
data$category <- factor(data$category, 
                        levels = c("Shared interaction relationships with ESCC SNPs", 
                                   "Unique interaction relationships"))

p1 <- ggplot(data, aes(x = "", y = value, fill = category)) +
  geom_bar(stat = "identity", width = 1) +
  coord_polar(theta = "y") +
  theme_void() +
  geom_text(aes(label = value), 
            position = position_stack(vjust = 0.5),
            size = 5) + 
  scale_fill_manual(values = c("Shared interaction relationships with ESCC SNPs" = "#D8BFD8", 
                               "Unique interaction relationships" = "#778899")) +
  ggtitle("Type 1") +
  theme(plot.title = element_text(hjust = 0.5, size = 14), 
        legend.position = "none",
        legend.text = element_text(size = 14))


#type2
data <- data.frame(
  category = c("Shared interaction relationships with ESCC SNPs", "Unique interaction relationships"),
  value = c(12, 93)
)
data$category <- factor(data$category, 
                        levels = c("Shared interaction relationships with ESCC SNPs", 
                                   "Unique interaction relationships"))

p2 <- ggplot(data, aes(x = "", y = value, fill = category)) +
  geom_bar(stat = "identity", width = 1) +
  coord_polar(theta = "y") +
  theme_void() +
  geom_text(aes(label = value), 
            position = position_stack(vjust = 0.5),
            size = 5) +  
  scale_fill_manual(values = c("Shared interaction relationships with ESCC SNPs" = "#D8BFD8", 
                               "Unique interaction relationships" = "#778899")) +
  ggtitle("Type 2") +
  theme(plot.title = element_text(hjust = 0.5, size = 14), 
        legend.position = "none",
        legend.text = element_text(size = 14))


#type3
data <- data.frame(
  category = c("Shared interaction relationships with ESCC SNPs", "Unique interaction relationships"),
  value = c(499, 2983)
)
data$category <- factor(data$category, 
                        levels = c("Shared interaction relationships with ESCC SNPs", 
                                   "Unique interaction relationships"))

p3 <- ggplot(data, aes(x = "", y = value, fill = category)) +
  geom_bar(stat = "identity", width = 1) +
  coord_polar(theta = "y") +
  theme_void() +
  geom_text(aes(label = value), 
            position = position_stack(vjust = 0.5),
            size = 5) + 
  scale_fill_manual(values = c("Shared interaction relationships with ESCC SNPs" = "#D8BFD8", 
                               "Unique interaction relationships" = "#778899"),
                    labels = c("Shared interaction\nrelationships\nwith ESCC SNPs",
                               "Unique interaction\nrelationships")) +
  ggtitle("Type 3") +
  theme(plot.title = element_text(hjust = 0.5, size = 14), 
        legend.title = element_blank(),
        legend.position = "right",         
        legend.key.height = unit(2, "cm"),
        legend.text = element_text(size = 14))

p <- p1 + p2 +p3
ggsave(filename = file.path("Fig3D-pie_2.png"), 
       plot = p, width = 10, height = 3, dpi = 300,bg = "white")
ggsave(filename = file.path("Fig3D-pie_2.pdf"), 
       plot = p, width = 10, height = 3, bg = "white")



#--------3E--------
#data preparation
#type1
cellline <- c("KYSE70","KYSE140","TT","KYSE180")
for(k in cellline){
  data <- read.table(paste0("type1_",k,".txt"))
  data$snpTarget <- paste(data$V4, data$V5, sep = "_"); data$cellLine <- k; data$value <- 1
  data <- data %>% select(snpTarget, cellLine, value)
  assign(paste0("type1_",k),data)
}
combined <- bind_rows(type1_KYSE70, type1_KYSE140, type1_TT, type1_KYSE180)
combined <- unique(combined)
result <- combined %>%
  select(snpTarget, cellLine) %>%   
  mutate(exists = 1) %>%                 
  pivot_wider(names_from = cellLine, 
              values_from = exists, 
              values_fill = 0) 
result$type <- "type1"
result$number <- result$KYSE140 + result$KYSE180 + result$KYSE70 + result$TT
result_type1 <- result %>% select(KYSE140,KYSE180, KYSE70, TT)
anno_type1 <- result %>% select(type, number)

#type2
cellline <- c("KYSE70","KYSE140","TT","KYSE180")
for(k in cellline){
  data <- read.table(paste0("type2_",k,".txt"))
  data$snpTarget <- paste(data$V4, data$V5, sep = "_"); data$cellLine <- k; data$value <- 1
  data <- data %>% select(snpTarget, cellLine, value)
  assign(paste0("type2_",k),data)
}
combined <- bind_rows(type2_KYSE70, type2_KYSE140, type2_TT, type2_KYSE180)
combined <- unique(combined)
result <- combined %>%
  select(snpTarget, cellLine) %>%   
  mutate(exists = 1) %>%                 
  pivot_wider(names_from = cellLine, 
              values_from = exists, 
              values_fill = 0) 
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

data <- rbind(result_type1,result_type2)
data <- t(data)
data <- as.matrix(data)
anno <- rbind(anno_type1, anno_type2)
anno <- t(anno)

top_anno = HeatmapAnnotation(type = anno[1,], number = anno[2,], border = T,
                             show_annotation_name = F,
                             col = list(
                               type = c("type1" = "darkcyan", "type2" = "slateblue4"),
                               number = c("1" = "turquoise1", "2" = "bisque", "3" = "darkorchid3","4" = "darkolivegreen1")
                             ),
                             annotation_legend_param = list(
                               type = list(title = "Types", 
                                           title_gp = gpar(fontsize = 14, fontface = "bold", col = "black"),
                                           labels_gp = gpar(fontsize = 14, col = "black")),
                               number = list(title = "Numbers", 
                                             title_gp = gpar(fontsize = 14, fontface = "bold", col = "black"),
                                             labels_gp = gpar(fontsize = 14, col = "black"))
                               
                             ))
p1 <- Heatmap(data,cluster_rows = F,
        cluster_columns = F,
        show_column_names = F,
        show_column_dend = F,
        col = c('white','deeppink2'),
        border = 'black',
        show_heatmap_legend = F,
        top_annotation =top_anno,
        row_names_side = "left",      # 将行名称显示在左侧
        row_names_gp = gpar(fontsize = 14))

data1 <- result_type3
data1 <- t(data1)
data1 <- as.matrix(data1)
anno1 <- anno_type3
anno1 <- t(anno1)

top_anno = HeatmapAnnotation(type = anno1[1,], number = anno1[2,], border = T,
                             show_annotation_name = F,
                             col = list(
                               type = c("type3" = "steelblue"),
                               number = c("1" = "turquoise1", "2" = "bisque", "3" = "darkorchid3","4" = "darkolivegreen1")
                             ),
                             annotation_legend_param = list(
                               type = list(title = "Types", 
                                           title_gp = gpar(fontsize = 14, fontface = "bold", col = "black"),
                                           labels_gp = gpar(fontsize = 14, col = "black")),
                               number = list(title = "Numbers", 
                                             title_gp = gpar(fontsize = 14, fontface = "bold", col = "black"),
                                             labels_gp = gpar(fontsize = 14, col = "black"))
                             ))
p2 <- Heatmap(data1,cluster_rows = F,
        cluster_columns = F,
        show_column_names = F,
        show_column_dend = F,
        col = c('white','deeppink2'),
        border = 'black',
        show_heatmap_legend = F,
        top_annotation =top_anno,
        row_names_side = "left",      # 将行名称显示在左侧
        row_names_gp = gpar(fontsize = 14))


g1 <- grid.grabExpr(draw(p1))
g2 <- grid.grabExpr(draw(p2))
p <- grid.arrange(g1, g2, ncol = 1)
p
ggsave(filename = file.path("Fig3E-snp_target.png"), 
       plot = p, width = 6.5, height = 4.5, dpi = 300,bg = "white")
ggsave(filename = file.path("Fig3E-snp_target.pdf"), 
       plot = p, width = 6.5, height = 4.5, bg = "white")



#--------3F--------
data <- data.frame(
  Category = c("type 1",  "type 2", "type 3", "type 3"),
  Group    = c("cis",  "cis", "cis", "trans"),
  Value    = c(45,  2,  75, 16)
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
ggsave(filename = file.path("Fig3F-eQTL_result.png"), 
       plot = plot1, width = 6, height = 5, dpi = 300,bg = "white")
ggsave(filename = file.path("Fig3F-eQTL_result.pdf"), 
       plot = plot1, width = 6, height = 5, bg = "white")
