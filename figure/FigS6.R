library(tidyverse)
library(pheatmap)

df <- read.csv(
  "TCN_incoming_delta_matrix.csv",
  check.names = FALSE
)


head(df)

heatmap_df <- df %>%
  column_to_rownames("LR_pair") %>%
  as.matrix()

heatmap_df <- apply(
  heatmap_df,
  2,
  as.numeric
)

rownames(heatmap_df) <- df$LR_pair

keep_lr <- apply(
  heatmap_df,
  1,
  function(x){
    sum(!is.na(x)) >= 3
  }
)


heatmap_df <- heatmap_df[keep_lr, ]

row_mean <- rowMeans(
  heatmap_df,
  na.rm = TRUE
)



heatmap_df <- heatmap_df[
  order(
    row_mean,
    decreasing = TRUE
  ),
]

max_value <- max(
  abs(heatmap_df),
  na.rm = TRUE
)


breaks <- seq(
  -0.005,
  0.005,
  length.out = 101
)


colors <- colorRampPalette(
  c(
    "#2166AC", 
    "#F7F7F7", 
    "#B2182B"   
  )
)(100)



pdf("heatmap_as_receiver.pdf",
    height = 4, width = 5.8)
pheatmap(
  heatmap_df,
  
  color = colors,
  
  breaks = breaks,
  
  na_col = "grey90",
  
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  
  border_color = NA,
  
  fontsize_row = 10,
  fontsize_col = 10,
  
  angle_col = 45,
  
  cellwidth = 25,
  cellheight = 18,
  
  main = "Differential LR communication score\n(TCN2-high vs TCN2-low TAM_LYVE1)"
)
dev.off()
