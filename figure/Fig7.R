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

#--------7A and 7B--------
Myeloid_sc_data <- readRDS("merge_seurat_Myeloid.RDS")
tam_cells <- subset(Myeloid_sc_data, subset = cell_type == "TAM_LYVE1")
tcn2_expr <- FetchData(tam_cells, vars = "TCN2")$TCN2

gmt_file <- "c2.cp.reactome.v2026.1.Hs.symbols.gmt"
pathways <- read.gmt(gmt_file)
pathway_list <- split(pathways$gene, pathways$term)
length(pathway_list)
min_genes <- 5
max_genes <- 500
pathway_list <- pathway_list[sapply(pathway_list, length) >= min_genes & 
                               sapply(pathway_list, length) <= max_genes]
cat("保留通路数量:", length(pathway_list), "\n")

tam_cells_ams <- tam_cells  
add_module_scores_safe <- function(obj, gene_sets, name_prefix = "AMS_", min_genes = 3) {
  score_list <- list()
  for (i in seq_along(gene_sets)) {
    pathway_name <- names(gene_sets)[i]
    genes <- gene_sets[[i]]
    available_genes <- intersect(genes, rownames(obj))
    if (length(available_genes) < min_genes) {
      warning("通路 ", pathway_name, " 只有 ", length(available_genes), 
              " 个基因存在于对象中，少于 ", min_genes, "，跳过")
      next
    }
    cat("处理 AddModuleScore 通路:", pathway_name, " (有效基因数:", length(available_genes), ")\n")
    obj <- AddModuleScore(obj, 
                          features = list(available_genes), 
                          name = paste0(name_prefix, pathway_name, "_"),
                          assay = "RNA")
    new_col <- paste0(name_prefix, pathway_name, "_1")
    if (new_col %in% colnames(obj@meta.data)) {
      score_list[[pathway_name]] <- obj@meta.data[[new_col]]
    } else {
      warning("通路 ", pathway_name, " 未成功添加评分（列未生成）")
    }
  }
  if (length(score_list) > 0) {
    scores_df <- as.data.frame(score_list)
    obj@meta.data <- cbind(obj@meta.data, scores_df)
  } else {
    scores_df <- data.frame()
    warning("没有成功添加任何通路评分")
  }
  return(list(obj = obj, scores = scores_df))
}

pathway_list_sub <- pathway_list   

res_ams <- add_module_scores_safe(tam_cells, pathway_list_sub, min_genes = 3)
tam_cells_ams <- res_ams$obj
ams_scores <- res_ams$scores

if (ncol(ams_scores) == 0) {
  stop("没有成功生成任何 AddModuleScore 评分，请检查通路列表和表达数据")
}

ams_cor <- data.frame(
  pathway = colnames(ams_scores),
  method = "AddModuleScore",
  rho = NA,
  p_value = NA,
  stringsAsFactors = FALSE
)

for (i in 1:ncol(ams_scores)) {
  ct <- cor.test(tcn2_expr, ams_scores[, i], method = "spearman")
  ams_cor$rho[i] <- ct$estimate
  ams_cor$p_value[i] <- ct$p.value
}


expr_mat <- GetAssayData(tam_cells, assay = "RNA", layer = "counts")
expr_mat <- as.matrix(expr_mat)  

cat("构建 AUCell 基因排序...\n")
cells_rankings <- AUCell_buildRankings(expr_mat, nCores = 4)


auc_scores_list <- list()
for (i in seq_along(pathway_list_sub)) {
  pathway_name <- names(pathway_list_sub)[i]
  cat("处理 AUCell 通路:", pathway_name, "\n")
  gene_set <- pathway_list_sub[[i]]
  gene_set <- intersect(gene_set, rownames(expr_mat))
  if (length(gene_set) < 3) {
    cat("通路", pathway_name, "在表达矩阵中基因过少，跳过\n")
    next
  }
  auc <- AUCell_calcAUC(setNames(list(gene_set), pathway_name), cells_rankings)
  auc_scores_list[[pathway_name]] <- as.numeric(getAUC(auc)[1, ])
}
auc_scores <- do.call(cbind, auc_scores_list)
colnames(auc_scores) <- names(auc_scores_list)

auc_cor <- data.frame(
  pathway = colnames(auc_scores),
  method = "AUCell",
  rho = NA,
  p_value = NA,
  stringsAsFactors = FALSE
)
for (i in 1:ncol(auc_scores)) {
  ct <- cor.test(tcn2_expr, auc_scores[, i], method = "spearman")
  auc_cor$rho[i] <- ct$estimate
  auc_cor$p_value[i] <- ct$p.value
}

ams_cor$p_adj <- p.adjust(ams_cor$p_value, method = "fdr")
auc_cor$p_adj <- p.adjust(auc_cor$p_value, method = "fdr")
sig_pathways <- intersect(
  ams_cor$pathway[ams_cor$p_adj < 0.05],
  auc_cor$pathway[auc_cor$p_adj < 0.05]
)

if (length(sig_pathways) > 0) {
  ams_sig <- ams_cor[ams_cor$pathway %in% sig_pathways, ]
  auc_sig <- auc_cor[auc_cor$pathway %in% sig_pathways, ]
  all_results <- rbind(ams_sig, auc_sig)
  all_results <- all_results[order(all_results$pathway, all_results$method), ]
  
  write.csv(all_results, "TCN2_pathway_correlation_significant.csv", row.names = FALSE)
  cat("分析完成！共有", length(sig_pathways), 
      "个通路在两种方法中均显著（FDR < 0.05）。结果已保存至 TCN2_pathway_correlation_significant.csv\n")
} else {
  cat("没有通路在两种方法中同时显著（FDR < 0.05）。\n")
}

sig_pathways <- all_results %>% filter(rho >= 0.40) %>% dplyr::select(pathway) %>% unique(.)
sig_pathways <- c(sig_pathways$pathway)

if (length(sig_pathways) > 0) {
  auc_scores_sig <- auc_scores[, colnames(auc_scores) %in% sig_pathways, drop = FALSE]
  all_genes <- rownames(expr_mat)
  
  all_cor_list <- list()
  
  for (pwy in colnames(auc_scores_sig)) {
    cat("正在处理通路:", pwy, " (计算所有基因与通路分数的相关性)...\n")
    pathway_score <- auc_scores_sig[, pwy]
    n_genes <- length(all_genes)
    res <- data.frame(
      pathway = rep(pwy, n_genes),
      gene = all_genes,
      rho = NA,
      p_value = NA,
      stringsAsFactors = FALSE
    )
    for (i in seq_along(all_genes)) {
      gene_expr <- expr_mat[all_genes[i], ]
      ct <- cor.test(gene_expr, pathway_score, method = "spearman")
      res$rho[i] <- ct$estimate
      res$p_value[i] <- ct$p.value
    }
    res$p_adj <- p.adjust(res$p_value, method = "fdr")
    res <- res[order(-abs(res$rho)), ]
    all_cor_list[[pwy]] <- res
  }
  
  all_genes_cor <- do.call(rbind, all_cor_list)
  write.csv(all_genes_cor, "all_genes_vs_significant_pathways.csv", row.names = FALSE)
  cat("所有基因与通路分数的相关性分析完成！结果保存至 all_genes_vs_significant_pathways.csv\n")
  top20_list <- list()
  for (pwy in unique(all_genes_cor$pathway)) {
    top20 <- head(all_genes_cor[all_genes_cor$pathway == pwy, ], 20)
    top20_list[[pwy]] <- top20
  }
  top20_all <- do.call(rbind, top20_list)
  write.csv(top20_all, "top20_genes_per_pathway.csv", row.names = FALSE)
  cat("每条通路 top 20 相关基因已保存至 top20_genes_per_pathway.csv\n")
  tcn2_res <- all_genes_cor[all_genes_cor$gene == "TCN2", ]
  if (nrow(tcn2_res) > 0) {
    cat("\nTCN2 在各显著通路中的相关性:\n")
    print(tcn2_res[, c("pathway", "rho", "p_adj")])
  }
  
} else {
  cat("没有显著通路，跳过全局相关性分析。\n")
}

target_pathway <- "REACTOME_TRANSPORT_OF_RCBL_WITHIN_THE_BODY"
pathway_score <- auc_scores[, target_pathway]
plot_df <- data.frame(
  TCN2_exp = tcn2_expr,
  Pathway_score = pathway_score
)
cor_test <- cor.test(plot_df$TCN2_exp, plot_df$Pathway_score, method = "spearman")
rho <- round(cor_test$estimate, 3)
p_val <- cor_test$p.value
p_label <- ifelse(p_val < 0.001, "p < 0.001", paste("p =", round(p_val, 3)))
p <- ggplot(plot_df, aes(x = TCN2_exp, y = Pathway_score)) +
  geom_point(alpha = 0.6, size = 1.5, color = "steelblue") +
  geom_smooth(method = "lm", se = TRUE, color = "red", linetype = "dashed",
              fill = "grey70", alpha = 0.3) +
  annotate("text", x = Inf, y = Inf, 
           label = paste("Spearman rho =", rho, "\n", p_label),
           hjust = 1.1, vjust = 1.1, size = 5, fontface = "plain") +  
  labs(x = "TCN2 Expression Level", 
       y = "AUCell Score:\nTRANSPORT OF RCBL WITHIN THE BODY") +
  theme_minimal(base_size = 12) +
  theme(
    panel.border = element_rect(fill = NA, color = "black", linewidth = 0.8),
    plot.title = element_text(hjust = 0.5),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.title = element_text(size = 14, color = "black"),  
    axis.text = element_text(size = 14, color = "black"),    
    legend.text = element_text(face = "plain")
  )
print(p)
ggsave(filename = file.path("Fig7A-1-TCN2_vs_B12_pathway_scatter.png"), 
       plot = p, width = 4, height = 4.5, dpi = 300,bg = "white")


target_pathway <- "REACTOME_ENDOSOMAL_VACUOLAR_PATHWAY"
pathway_score <- auc_scores[, target_pathway]
plot_df <- data.frame(
  TCN2_exp = tcn2_expr,
  Pathway_score = pathway_score
)
cor_test <- cor.test(plot_df$TCN2_exp, plot_df$Pathway_score, method = "spearman")
rho <- round(cor_test$estimate, 3)
p_val <- cor_test$p.value
p_label <- ifelse(p_val < 0.001, "p < 0.001", paste("p =", round(p_val, 3)))
p <- ggplot(plot_df, aes(x = TCN2_exp, y = Pathway_score)) +
  geom_point(alpha = 0.6, size = 1.5, color = "steelblue") +
  geom_smooth(method = "lm", se = TRUE, color = "red", linetype = "dashed",
              fill = "grey70", alpha = 0.3) +
  annotate("text", x = Inf, y = Inf, 
           label = paste("Spearman rho =", rho, "\n", p_label),
           hjust = 1.1, vjust = 1.1, size = 5, fontface = "plain") + 
  labs(x = "TCN2 Expression Level", 
       y = "AUCell Score:\nENDOSOMAL VACUOLAR PATHWAY") +
  theme_minimal(base_size = 12) +
  theme(
    panel.border = element_rect(fill = NA, color = "black", linewidth = 0.8),
    plot.title = element_text(hjust = 0.5),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 12),
    legend.text = element_text(face = "plain")
  )
print(p)
ggsave(filename = file.path("Fig7B-2-TCN2_vs_ENDOSOMAL_pathway_scatter.png"), 
       plot = p, width = 4, height = 4.5, dpi = 300,bg = "white")


HLAE_expr <- FetchData(tam_cells, vars = "HLA-E")$`HLA-E`
target_pathway <- "REACTOME_ENDOSOMAL_VACUOLAR_PATHWAY"
pathway_score <- auc_scores[, target_pathway]
plot_df <- data.frame(
  HLAE_expr = HLAE_expr,
  Pathway_score = pathway_score
)
cor_test <- cor.test(plot_df$HLAE_expr, plot_df$Pathway_score, method = "spearman")
rho <- round(cor_test$estimate, 3)
p_val <- cor_test$p.value
p_label <- ifelse(p_val < 0.001, "p < 0.001", paste("p =", round(p_val, 3)))
p <- ggplot(plot_df, aes(x = HLAE_expr, y = Pathway_score)) +
  geom_point(alpha = 0.6, size = 1.5, color = "steelblue") +
  geom_smooth(method = "lm", se = TRUE, color = "red", linetype = "dashed",
              fill = "grey70", alpha = 0.3) +
  annotate("text", x = Inf, y = Inf, 
           label = paste("Spearman rho =", rho, "\n", p_label),
           hjust = 1.1, vjust = 1.1, size = 5, fontface = "plain") + 
  labs(x = "HLA-E Expression Level", 
       y = "AUCell Score:\nENDOSOMAL VACUOLAR PATHWAY") +
  theme_minimal(base_size = 12) +
  theme(
    panel.border = element_rect(fill = NA, color = "black", linewidth = 0.8),
    plot.title = element_text(hjust = 0.5),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 12),
    legend.text = element_text(face = "plain")
  )
print(p)
ggsave(filename = file.path("Fig7B-3-HLAE_vs_B12_pathway_scatter.png"), 
       plot = p, width = 4, height = 4.5, dpi = 300,bg = "white")

#--------7C--------
output <- file.path("")
load(file = file.path(output, "cellchat_ASS-ILD.RData"))
png(file.path("Fig7C-cellchat.png"), width = 1500, height = 1500, res = 300)
netVisual_chord_gene(cellchat_AS, sources.use = 16, legend.pos.x = 2)
dev.off()

pdf(file.path("Fig7C-cellchat.pdf"), width = 5, height = 5)
netVisual_chord_gene(cellchat_AS, sources.use = 16, legend.pos.x = 2)
dev.off()


#--------7D and 7E--------
output <- file.path("")
load(file = file.path(output, "mycds.Rdata"))

p1 <- plot_cell_trajectory(mycds, color_by = "Pseudotime", show_backbone = TRUE, cell_size = 0.5) +
  theme(
    legend.position = "right",                     
    axis.title = element_text(size = 15.5, color = "black"),
    axis.text = element_text(size = 14, color = "black"),
    legend.text = element_text(size = 15.5),
    legend.title = element_text(size = 15.5)
  ) +
  scale_color_viridis_c() +
  guides(
    color = guide_colorbar(
      barwidth = 1,                               
      barheight = 10,                         
      title.position = "top",
      label.theme = element_text(size = 15.5)
    )
  )
p1

p4 <- plot_cell_trajectory(mycds, color_by = "cell_type", show_backbone = TRUE, cell_size = 0.5) +
  theme(
    legend.position = "right",  
    axis.title = element_text(size = 15.5, color = "black"),
    axis.text = element_text(size = 14, color = "black"),
    legend.text = element_text(size = 15.5),
    legend.title = element_text(size = 15.5)
  ) +
  guides(
    color = guide_legend(
      ncol = 1,               
      byrow = TRUE,
      override.aes = list(size = 3)
    )
  )
p4
combined_plot <- p1 + p4
combined_plot
ggsave(file.path("Fig7D_Trajectory_Plot.png"), combined_plot,
       width = 9, height = 4, dpi = 300, bg = "white")
ggsave(file.path("Fig7D_Trajectory_Plot.pdf"), combined_plot,
       width = 9, height = 4,bg = "white")
p <- plot_genes_in_pseudotime(mycds[c("KLRC1"),],
                              color_by = "cell_type",
                              nrow = 1, 
                              ncol = NULL) +
  theme(
    legend.position = "none",  
    axis.title = element_text(size = 15.5, color = "black"),
    axis.text = element_text(size = 15.5, color = "black"),
    strip.text = element_text(size = 15.5, color = "black")  
  )
p
ggsave(file.path("Fig7E_gene_expresssion.png"), p,
       width = 3, height = 4, dpi = 300, bg = "white")
ggsave(file.path("Fig7E_gene_expresssion.pdf"), p,
       width = 3, height = 4, bg = "white")


