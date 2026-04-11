library(devtools)
library(ComplexHeatmap)
library(ggplot2)
library(trackViewer)
library(dplyr)
library(tidyr)

rm(list = ls()); gc()
ORIGINAL_DIR <- ""
setwd(ORIGINAL_DIR)

#--------1B--------
#chr2
SNP <- c(135677161, 192207115, 202122995, 202162811, 230919172)
SNPname <- c(12, 13, 14, 15, 16)
sample.gr <- GRanges("chr2", IRanges(SNP, width=1, names=paste0("snp", SNPname)))
sample.gr$color <- sample.int(6, length(SNP), replace=TRUE)

snp_colors <- c("red", "green", "purple", "orange","yellow") 
sample.gr$color <- snp_colors

genename <- c("CCNT2","MYO1B","CASP8","ALS2CR12","SLC16A14")
features <- GRanges("chr2", IRanges(c(135675804, 192109910, 202098165, 202152993,230899697), 
                                    width=c(41108,180205,54269,69128,34018),
                                    names=genename))
features$fill <- c("#FF8833", "#51C6E6", "#DFA32D", "#FF7893", "#496")
features$height <- c(0.08, 0.08, 0.08, 0.08, 0.08)

chr2 <- sample.gr
chr2$label.parameter.rot <- 45

xaxis <- c(135675804, 192109910, 202098165, 202152993,230933715)
rescale <- data.frame(
  from.start = c(135675804,135716913,192109910,192290116,202098165,202152435,202152993,202222121,230899697), 
  from.end   = c(135716912,192109909,192290115,202098164,202152434,202152992,202222121,230899696,230933715),
  to.start   = c(135675804,143947506,155271126,191531071,193500511,204420387,204420499,218330212,224088622), 
  to.end     = c(143947505,155271125,191531070,193500510,204420386,204420498,218330211,224088621,230933715)
)
png(file.path("Fig1B-location_chr2.png"), width = 2000, height = 600, res = 300)
lolliplot(chr2, features,yaxis=FALSE,type="pin",cex=1.3,xaxis=xaxis, rescale=rescale,
          label.parameter.color=chr2$color,label.parameter.cex = 3,
          legend.parameter.offset = 0)
dev.off()
pdf(file.path("Fig1B-location_chr2.pdf"), width = 7, height = 2)
lolliplot(chr2, features,yaxis=FALSE,type="pin",cex=1.3,xaxis=xaxis, rescale=rescale,
          label.parameter.color=chr2$color,label.parameter.cex = 3,
          legend.parameter.offset = 0)
dev.off()

#chr9
SNP <- c(21974218, 21997015, 22033366)
SNPname <- c(28,29,30)
sample.gr <- GRanges("chr9", IRanges(SNP, width=1, names=paste0("snp", SNPname)))
sample.gr$color <- sample.int(6, length(SNP), replace=TRUE)

snp_colors <- c("red", "green", "orange")
sample.gr$color <- snp_colors

genename <- c("CDKN2A","CDKN2B-AS1","UBA52P6")
features <- GRanges("chr9", IRanges(c(21967750, 21994776,22012535), 
                                    width=c(27750,126320,382),
                                    names=genename))
features$fill <- c("#FF8833", "#51C6E6", "#FF7893")
features$height <- c(0.10, 0.08, 0.10)

chr9 <- sample.gr
chr9$label.parameter.rot <- 45

xaxis <- c(21967750,22012535,22121096)
rescale <- data.frame(
  from.start = c(21967750,21994776,21995301,22012154,22012536), 
  from.end   = c(21994775,21995300,22012153,22012535,22121096),
  to.start   = c(21967750,21994776,21995301,22012154,22012536), 
  to.end     = c(21994775,21995300,22012153,22012535,22121096)
)

png(file.path("Fig1B-location_chr9.png"), width = 2000, height = 600, res = 300)
lolliplot(chr9, features,yaxis=FALSE,type="pin",cex=1.3,xaxis=xaxis, rescale=rescale,
          label.parameter.color=chr9$color,label.parameter.cex = 3,
          legend.parameter.offset = 0)
dev.off()
pdf(file.path("Fig1B-location_chr9.pdf"), width = 7, height = 2)
lolliplot(chr9, features,yaxis=FALSE,type="pin",cex=1.3,xaxis=xaxis, rescale=rescale,
          label.parameter.color=chr9$color,label.parameter.cex = 3,
          legend.parameter.offset = 0)
dev.off()

#chr22
SNP <- c(29156448,29198151,29375026,30986350)
SNPname <- c(17,18,19,20)
sample.gr <- GRanges("chr22", IRanges(SNP, width=1, names=paste0("snp", SNPname)))
sample.gr$color <- sample.int(6, length(SNP), replace=TRUE)

snp_colors <- c("red", "green", "yellow", "orange")
sample.gr$color <- snp_colors

genename <- c("HSCB","CCDC117","XBP1","ZNRF3","PES1")
features <- GRanges("chr22", IRanges(c(29138018,29168661,29190542,29279579,30972611), 
                                     width=c(15485,16622,6043,173896,30459),
                                     names=genename))
features$fill <- c("#FF8833", "#51C6E6", "#DFA32D", "#FF7893", "#496")
features$height <- c(0.08, 0.08, 0.08, 0.08, 0.08)

chr22 <- sample.gr
chr22$label.parameter.rot <- 45

xaxis <- c(29138018,29168661,29190542,29279579,31003070)
rescale <- data.frame(
  from.start = c(29138018,29153504,29168661,29185284,29190542,29196586,29279579,29453476,30972611), 
  from.end   = c(29153503,29168660,29185283,29190541,29196585,29279578,29453475,30972610,31003070),
  to.start   = c(29138018,29233295,29236780,29339053,29340262,29377448,29396527,30466428,30815666), 
  to.end     = c(29233294,29236779,29339052,29340261,29377447,29396526,30466427,30815665,31003070)
)
png(file.path("Fig1B-location_chr22.png"), width = 2000, height = 600, res = 300)
lolliplot(chr22, features,yaxis=FALSE,type="pin",cex=1.3,xaxis=xaxis, rescale=rescale,
          label.parameter.color=chr22$color,label.parameter.cex = 3,
          legend.parameter.offset = 0)
dev.off()

pdf(file.path("Fig1B-location_chr22.pdf"), width = 7, height = 2)
lolliplot(chr22, features,yaxis=FALSE,type="pin",cex=1.3,xaxis=xaxis, rescale=rescale,
          label.parameter.color=chr22$color,label.parameter.cex = 3,
          legend.parameter.offset = 0)
dev.off()

#--------1D--------
snp <-read.csv(file='SNP.CSV',
               header=TRUE,sep=',')
snp<- as.matrix(snp)

anno <- read.csv(file='SNP1.CSV',
                 header=TRUE,sep=',')
left_anno=rowAnnotation(df=anno,border = T,
                        show_annotation_name = F,
                        col = list(
                          Types = c("promoter" = "grey", "peaks" = "skyblue", "anchor peaks" = "tan1"),
                          Factors = c("H3K27ac"="goldenrod1", "H3K4me1"="slateblue1", "H3K4me3"="pink", "promoter"="palegreen3", "SOX2"="lightseagreen"),
                          CellLines = c("all"="darkcyan", "KYSE140"="chocolate3","KYSE180"="black","KYSE70"="wheat3","TE5"="purple3","TT"="blue")
                        ),
                        annotation_legend_param = list(
                          Types = list(title = "Types", 
                                      title_gp = gpar(fontsize = 14, fontface = "bold", col = "black"),
                                      labels_gp = gpar(fontsize = 14, col = "black")),
                          Factors = list(title = "Factors", 
                                      title_gp = gpar(fontsize = 14, fontface = "bold", col = "black"),
                                      labels_gp = gpar(fontsize = 14, col = "black")),
                          CellLines = list(title = "Cell Line", 
                                      title_gp = gpar(fontsize = 14, fontface = "bold", col = "black"),
                                      labels_gp = gpar(fontsize = 14, col = "black"))
                        ))

png(file.path("Fig1D-heatmap_snp.png"), width = 2300, height = 1400, res = 300)
Heatmap(snp,cluster_rows = F,
        cluster_columns = T,
        show_column_names = F,
        show_column_dend = FALSE,
        col = c('white','black'),
        border = 'black',
        show_heatmap_legend = F,
        left_annotation =left_anno )
dev.off()

pdf(file.path("Fig1D-heatmap_snp.pdf"), width = 7.5, height = 5)
Heatmap(snp,cluster_rows = F,
        cluster_columns = T,
        show_column_names = F,
        show_column_dend = FALSE,
        col = c('white','black'),
        border = 'black',
        show_heatmap_legend = F,
        left_annotation =left_anno )
dev.off()


#--------1F--------
data2 <- read.csv(file='heat2.CSV',
                  header=TRUE,sep=',')
data3 <- read.csv(file='heat3.csv',
                  header=TRUE,sep=',')

p1 <- ggplot(data2, aes(cellline,factor))+
  geom_point(aes(size = effectsize,fill=P),shape=21,stroke=0.6) +
  scale_fill_gradient(name = 'P',
                      limit = c(-0.00001,0.050001),
                      breaks = c(0.001,0.01,0.05),
                      low = 'deeppink',
                      high = 'white')+
  scale_size_continuous(name = 'Effect Size',
                        limit = c(-0.001,50),
                        breaks = c(0,5,10,15),
                        range = c(3,17))+
  
  scale_y_discrete(limits=as.character(c("SOX2","H3K4me3","H3K4me1","H3K27ac","H3K27ac anchors")))+ 
  scale_x_discrete(limits=as.character(c("KYSE140","KYSE180","KYSE70","TT","TE5")))+ 
  theme_bw()+
  xlab(NULL) + 
  ylab(NULL)+
  theme(panel.border = element_rect(fill=NA,color="black", size=0.5, linetype="solid"),
        axis.text.x = element_text(size = 14, angle = 45, hjust = 1,colour = "black"),  
        axis.text.y = element_text(size = 14, angle = 45,colour = "black"),  
        legend.title = element_text(size = 14),  
        legend.text = element_text(size = 12))  
p1
p2<-p1+geom_point(data=data3,
                  mapping =aes(cellline,factor),
                  shape=4,
                  stroke=1,
                  size=8,
                  color='black')
p2
ggsave(filename = file.path("Fig1F-shuffle.png"), 
       plot = p2, width = 5, height = 4.3, dpi = 300,bg = "white")
ggsave(filename = file.path("Fig1F-shuffle.pdf"), 
       plot = p2, width = 5, height = 4.3, bg = "white")

#--------1G--------
data2 <- read.csv(file='heat2SNP.CSV',
                  header=TRUE,sep=',')
data3 <- read.csv(file='heat3.csv',
                  header=TRUE,sep=',')

p1 <- ggplot(data2, aes(cellline,factor))+
  geom_point(aes(size = effectsize,fill=P),shape=21,stroke=0.6) +
  scale_fill_gradient(name = 'P',
                      limit = c(-0.00001,0.050001),
                      breaks = c(0.001,0.01,0.05),
                      low = 'deeppink',
                      high = 'white')+
  scale_size_continuous(name = 'Effect Size',
                        limit = c(-0.001,50),
                        breaks = c(0,5,10,15),
                        range = c(3,17))+
  
  scale_y_discrete(limits=as.character(c("SOX2","H3K4me3","H3K4me1","H3K27ac","H3K27ac anchors")))+ 
  scale_x_discrete(limits=as.character(c("KYSE140","KYSE180","KYSE70","TT","TE5")))+ 
  theme_bw()+
  xlab(NULL) + 
  ylab(NULL)+
  theme(panel.border = element_rect(fill=NA,color="black", size=0.5, linetype="solid"),
        axis.text.x = element_text(size = 14, angle = 45, hjust = 1,colour = "black"),   
        axis.text.y = element_text(size = 14, angle = 45,colour = "black"),   
        legend.title = element_text(size = 14),  
        legend.text = element_text(size = 12))  
p1
p2<-p1+geom_point(data=data3,
                  mapping =aes(cellline,factor),
                  shape=4,
                  stroke=1,
                  size=8,
                  color='black')
p2
ggsave(filename = file.path("Fig1G-shuffle_snp.png"), 
       plot = p2, width = 5, height = 4.3, dpi = 300,bg = "white")
ggsave(filename = file.path("Fig1G-shuffle_snp.pdf"), 
       plot = p2, width = 5, height = 4.3, bg = "white")
