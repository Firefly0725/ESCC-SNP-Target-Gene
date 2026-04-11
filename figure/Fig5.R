library(limma)
library(ggplot2)
library(pheatmap)
library(ggrepel)
library(ComplexHeatmap)
library(circlize)
library(survival)
library(ggtext)
library(survminer)

rm(list = ls()); gc()
ORIGINAL_DIR <- ""
setwd(ORIGINAL_DIR)

#--------5A--------
data_dir <- "E:/ESCC_SNP"
gse_file <- file.path(data_dir, "GSE53625_clean.RData")

de_adjP_cutoff <- 0.05
de_logFC_cutoff <- 0.5
plot_base_size <- 14

loaded_names <- load(gse_file)
gse <- GSE53625_clean

expr_raw <- as.matrix(gse$expr)          
clinical <- as.data.frame(gse$clinical, stringsAsFactors = FALSE)

common <- intersect(colnames(expr_raw), rownames(clinical))
expr <- expr_raw[, common, drop = FALSE]
clinical <- clinical[common, , drop = FALSE]

group_raw <- trimws(tolower(clinical$group))
group <- ifelse(group_raw %in% c("tumor","tumour","cancer","escc","case"), "cancer",
                ifelse(group_raw %in% c("normal","control","adjacent normal","adjacent_normal"), "normal", NA))

n <- ncol(expr)
cancer_pos <- seq(1, n, by = 2)
normal_pos <- seq(2, n, by = 2)

pair_id <- factor(rep(1:(n/2), each = 2))
tissue <- factor(group, levels = c("normal", "cancer"))

pair_info <- data.frame(sample_id = colnames(expr), pair_id = pair_id, tissue = tissue)

design <- model.matrix(~ pair_id + tissue, data = pair_info)
fit <- lmFit(expr, design)
fit <- eBayes(fit)

coef_name <- "tissuecancer"
de_all <- topTable(fit, coef = coef_name, number = Inf, sort.by = "P", adjust.method = "BH")
de_all$gene <- rownames(de_all)
rownames(de_all) <- NULL

de_sig <- de_all[de_all$adj.P.Val < de_adjP_cutoff & abs(de_all$logFC) >= de_logFC_cutoff, ]


de_all$neglog10FDR <- -log10(pmax(de_all$adj.P.Val, 1e-300))
de_all$DE_status <- "NS"
de_all$DE_status[de_all$adj.P.Val < de_adjP_cutoff & de_all$logFC >  de_logFC_cutoff] <- "Up"
de_all$DE_status[de_all$adj.P.Val < de_adjP_cutoff & de_all$logFC < -de_logFC_cutoff] <- "Down"


load("filtered_biomarker.Rdata")
genes_to_label <- biomarker
de_all_sig <- de_all %>% filter(DE_status %in% c("Up","Down"))

label_data <- de_all_sig[de_all_sig$gene %in% genes_to_label, ]

volcano_plot <- ggplot(de_all, aes(x = logFC, y = neglog10FDR)) +
  geom_point(aes(color = DE_status), size = 1.2, alpha = 0.75) +
  geom_vline(xintercept = c(-de_logFC_cutoff, de_logFC_cutoff), linetype = 2, color = "grey60") +
  geom_hline(yintercept = -log10(de_adjP_cutoff), linetype = 2, color = "grey60") +
  scale_color_manual(values = c(NS = "grey70", Up = "#D55E00", Down = "#0072B2")) +
  geom_point(data = label_data, 
             aes(x = logFC, y = neglog10FDR),
             color = "black", fill = "gold", shape = 21, size = 3, stroke = 0.8) +
  geom_text_repel(data = label_data,
                  aes(x = logFC, y = neglog10FDR, label = gene),
                  size = 14 / .pt,   
                  box.padding = 0.35, point.padding = 0.2,
                  segment.color = "black", segment.size = 0.5) +
  labs(x = "logFC (Tumor - Normal)", y = expression(-log[10]("FDR")), color = NULL) +
  theme_classic(base_size = plot_base_size) +  
  theme(plot.title = element_blank(), 
        panel.grid = element_blank(),
        axis.text = element_text(color = "black", size = 14),   
        axis.title = element_text(color = "black", size = 14))
volcano_plot
ggsave(file.path("Fig5A-paired_DE_volcano.png"), volcano_plot, width = 5.5, height = 5, dpi = 300)
ggsave(file.path("Fig5A-paired_DE_volcano.pdf"), volcano_plot, width = 5.5, height = 5, bg = "white")


#--------S1A--------
font_size <- 14
row_var <- apply(expr, 1, var, na.rm = TRUE)
top_genes <- names(sort(row_var, decreasing = TRUE))[1:min(30, nrow(expr))]
mat <- expr[top_genes, , drop = FALSE]
mat_z <- t(scale(t(mat)))
mat_z[!is.finite(mat_z)] <- 0

tissue_colors <- c("cancer" = "#D55E00", "normal" = "#0072B2")
pair_levels <- levels(pair_info$pair_id)
n_pair <- length(pair_levels)
pair_colors <- setNames(rainbow(n_pair), pair_levels)  # 或使用 viridis

annotation_col <- data.frame(
  Tissue = pair_info$tissue,
  Pair = pair_info$pair_id
)
rownames(annotation_col) <- pair_info$sample_id

ha <- HeatmapAnnotation(
  df = annotation_col,
  col = list(
    Tissue = tissue_colors,
    Pair = pair_colors
  ),
  show_legend = c(Tissue = FALSE, Pair = FALSE)
)

all_genes <- rownames(mat_z)
row_labels <- ifelse(all_genes %in% label_data$gene, all_genes, "")

ht <- Heatmap(
  mat_z,
  name = "Z-score",
  col = colorRamp2(c(-2, 0, 2), c("blue", "white", "red")),
  top_annotation = ha,
  show_column_names = FALSE,
  row_labels = row_labels,
  row_names_gp = gpar(fontsize = font_size),
  column_title_gp = gpar(fontsize = font_size),
  heatmap_legend_param = list(
    title = "Z-score",
    at = c(-2, 0, 2),
    title_gp = gpar(fontsize = font_size),
    labels_gp = gpar(fontsize = font_size)
  )
)

tissue_legend <- Legend(
  title = "Tissue",
  at = c("cancer", "normal"),
  legend_gp = gpar(fill = tissue_colors),
  labels = c("Cancer", "Normal"),
  title_gp = gpar(fontsize = font_size),
  labels_gp = gpar(fontsize = font_size)
)

draw(ht,
     annotation_legend_list = list(tissue_legend),
     annotation_legend_side = "right",
     merge_legend = FALSE)


png("SFig1A-heatmap_with_tissue_legend.png", width = 2000, height = 1600, res = 300)
draw(ht,
     annotation_legend_list = list(tissue_legend),
     annotation_legend_side = "right",
     merge_legend = FALSE)
dev.off()

pdf("SFig1A-heatmap_with_tissue_legend.pdf", width = 7, height = 5)
draw(ht,
     annotation_legend_list = list(tissue_legend),
     annotation_legend_side = "right",
     merge_legend = FALSE)
dev.off()



#--------5B--------
data_dir <- "E:/ESCC_SNP"
discovery_file   <- file.path("GSE53625_clean_cancer.RData")
replication_file <- file.path("ESCC_two_tcga_extracted_clean.RData")
replication_object_name <- "ESCC_tcga_legacy"
biomarker_file <- file.path(data_dir, "filtered_gene.Rdata")

discovery_p_cutoff <- 0.05
use_discovery_fdr  <- TRUE
discovery_fdr_cutoff <- 0.10
require_ph_ok <- TRUE
ph_p_cutoff <- 0.05
forest_top_n <- 25
plot_base_size <- 15

load(file = "Fig5_dis_and_val_expr_clin.Rdata")
load(file = "discovery_univariate_cox_all.Rdata")

plot_df <- head(discovery_res, forest_top_n)
plot_df$gene <- factor(plot_df$gene, levels = rev(plot_df$gene))
  
theme_pub <- theme_classic(base_size = plot_base_size) +
    theme(plot.title = element_blank(),
          panel.grid = element_blank(),
          axis.title = element_text(size = plot_base_size, face = "plain"),
          axis.text  = element_text(size = plot_base_size, face = "plain"),
          legend.title = element_text(size = plot_base_size, face = "plain"),
          legend.text  = element_text(size = plot_base_size, face = "plain"))
  
p <- ggplot(plot_df, aes(x = HR, y = gene)) +
  geom_vline(xintercept = 1, linetype = 2, color = "grey60", linewidth = 0.8) +     
  geom_errorbarh(aes(xmin = HR_L, xmax = HR_H), height = 0.20, color = "grey9", linewidth = 0.8) +  
  geom_point(size = 3, color = "#D55E00") +
  scale_x_log10() +
  labs(x = "Hazard ratio (per 1 SD increase, log scale)", y = NULL) +
  scale_y_discrete(labels = function(x) ifelse(x == "TCN2", paste0("**", x, "**"), x)) +
  theme_pub +
  theme(axis.text.y = ggtext::element_markdown(size = plot_base_size, face = "plain", color = "black"),
        axis.text.x = element_text(size = plot_base_size, face = "plain", color = "black"))
p 
ggsave(file.path("Fig5B-discovery_top_forestplot.png"), p, width = 7,
         height = 8, dpi = 300)
ggsave(file.path("Fig5B-discovery_top_forestplot.pdf"), p, width = 7,
       height = 8, bg = "white")



#--------5C--------
plot_base_size = 15
load(file = "Fig5C-rep_df.Rdata")
long_df <- rbind(
  data.frame(
    gene = rep_df$gene,
    cohort = "Discovery",
    HR = rep_df$HR_discovery,
    HR_L = rep_df$HR_L_discovery,
    HR_H = rep_df$HR_H_discovery
  ),
  data.frame(
    gene = rep_df$gene,
    cohort = "Replication",
    HR = rep_df$HR_replication,
    HR_L = rep_df$HR_L_replication,
    HR_H = rep_df$HR_H_replication
  )
)
long_df$gene <- factor(long_df$gene, levels = rev(unique(rep_df$gene)))
long_df$cohort <- factor(long_df$cohort, levels = c("Discovery", "Replication"))
  
theme_pub <- theme_classic(base_size = plot_base_size) +
  theme(plot.title = element_blank(),
          panel.grid = element_blank(),
          axis.title = element_text(size = plot_base_size, face = "plain"),
          axis.text = element_text(size = plot_base_size, face = "plain"),
          legend.title = element_text(size = plot_base_size, face = "plain"),
          legend.text = element_text(size = plot_base_size, face = "plain"))
p <- ggplot(long_df, aes(x = HR, y = gene, color = cohort)) +
  geom_vline(xintercept = 1, linetype = 2, color = "grey60", linewidth = 1) +
  geom_errorbarh(aes(xmin = HR_L, xmax = HR_H), height = 0.20, 
                 position = position_dodge(width = 0.55), linewidth = 1) +
  geom_point(position = position_dodge(width = 0.55), size = 2.5) +
  scale_x_log10() +
  scale_color_manual(values = c("Discovery" = "#D55E00",   
                                "Replication" = "#0072B2"), 
                     name = "Cohort") +
  labs(x = "Hazard ratio (per 1 SD increase, log scale)", y = NULL, color = NULL) +
  theme_pub +
  theme(axis.text = element_text(color = "black")) 
p
ggsave(file.path("Fig3C-replicated_gene_forestplot.png"), p, width = 5.5,
       height = 2.5, dpi = 300)
ggsave(file.path("Fig3C-replicated_gene_forestplot.pdf"), p, width = 5.5,
       height = 2.5)

#--------5D and 5E--------
data_dir <- "E:/ESCC_SNP"

discovery_file  <- file.path("GSE53625_clean_cancer.RData")
validation_file <- file.path("ESCC_two_tcga_extracted_clean.RData")
discovery_label  <- "Discovery: GSE53625"
validation_label <- "Validation: TCGA legacy ESCC"
load(file = "Fig5C-rep_df.Rdata")
merged_results_file <- rep_df
gene_to_plot <- "TCN2"   
split_method <- "median"
min_group_size <- 5
cex_base <- 1.3


merged_df <- merged_results_file
if (!is.null(gene_to_plot) && gene_to_plot %in% merged_df$gene) {
  gene_selected <- gene_to_plot
} else if (any(merged_df$replicated, na.rm = TRUE)) {
  tmp <- merged_df[merged_df$replicated, ]
  gene_selected <- tmp$gene[order(tmp$replication_pval, tmp$discovery_pval)][1]
} else if (any(merged_df$direction_consistent, na.rm = TRUE)) {
  tmp <- merged_df[merged_df$direction_consistent, ]
  gene_selected <- tmp$gene[order(tmp$replication_pval, tmp$discovery_pval)][1]
} else {
  gene_selected <- merged_df$gene[order(merged_df$discovery_pval, merged_df$replication_pval)][1]
}
message("Selected gene: ", gene_selected)

load( file = "Fig5_dis_and_val_expr_clin.Rdata")

plot_km_with_risk <- function(expr_mat, clin_df, gene, cohort_label) {
  if (!gene %in% rownames(expr_mat)) return(NULL)
  x <- as.numeric(expr_mat[gene, ])
  keep <- is.finite(x) & is.finite(clin_df$OS.time) & is.finite(clin_df$OS.event)
  x <- x[keep]; time <- clin_df$OS.time[keep]; event <- clin_df$OS.event[keep]
  if (length(unique(x)) < 2) return(NULL)
  cutv <- median(x, na.rm = TRUE)
  group <- factor(ifelse(x > cutv, "High", "Low"), levels = c("Low", "High"))
  if (any(table(group) < 5)) return(NULL)
  fit <- survfit(Surv(time, event) ~ group)
  
  my_theme <- theme_bw() +
    theme(panel.grid = element_blank(),
          axis.title = element_text(size = 15),
          axis.text = element_text(size = 15),
          legend.title = element_text(size = 15),
          legend.text = element_text(size = 15),
          plot.title = element_text(size = 15))
  my_table_theme <- theme_cleantable() +
    theme(axis.text = element_text(size = 15),
          text = element_text(size = 15))
  p <- ggsurvplot(fit, 
                  data = data.frame(time = time, event = event, group = group),
                  risk.table = TRUE,
                  risk.table.col = "strata",
                  pval = TRUE,
                  pval.method = FALSE,
                  conf.int = FALSE,
                  palette = c("#0072B2", "#D55E00"),
                  xlab = "Overall survival time",
                  ylab = "Survival probability",
                  legend.title = cohort_label,
                  legend.labs = c("Low", "High"),
                  risk.table.y.text.col = TRUE,
                  risk.table.y.text = FALSE,
                  tables.theme = my_table_theme,
                  ggtheme = my_theme,
                  font.x = 15,
                  font.y = 15,
                  font.tickslab = 15,
                  font.legend = 15,
                  font.pval = 15,
                  font.risk.table = 15)
  if (!is.null(p$table)) {
    time_range <- range(c(0, time), na.rm = TRUE)
    breaks <- pretty(time_range)
    p$table <- p$table +
      xlab("Time") +
      scale_x_continuous(breaks = breaks, limits = time_range) +
      theme(axis.text.x = element_text(size = 15),
            axis.title.x = element_text(size = 15))
  }
  return(p)
}

png(file.path("Fig5D-KM-dis.png"), width=1600, height=1600, res=300)
plot_km_with_risk(expr_dis, clin_dis, gene_selected, discovery_label)
dev.off()

png(file.path("Fig5E-KM-val.png"), width=1600, height=1600, res=300)
plot_km_with_risk(expr_val, clin_val, gene_selected, validation_label)
dev.off()

pdf(file.path("Fig5D-KM-dis.pdf"), width=6, height=6)
plot_km_with_risk(expr_dis, clin_dis, gene_selected, discovery_label)
dev.off()

pdf(file.path("Fig5E-KM-val.pdf"), width=6, height=6)
plot_km_with_risk(expr_val, clin_val, gene_selected, validation_label)
dev.off()


