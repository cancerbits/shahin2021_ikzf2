## ---- parameter_list
blueprint.se <- BlueprintEncodeData()
workers <- 16
BPPARAM <- MulticoreParam(workers=workers)
theme_set(theme_cowplot()) #ggplot2 theme

seuratAnalysisCfg <- list(
  PCA_DIMS = 1:30,
  MIN_FEATURES = 500,
  MIN_CELLS = 20,
  PERC_MITO = 10,
  MAX_DOUBLET_SCORE = 3
  )

diff_test <- "MAST"
min_cell_n <- 40 
lfc_threshold <- 0.5
fdr_threshold <- 0.05

## ---- colors
pt_cols <- c(brewer.pal(4,"Dark2"),
             brewer.pal(8,"Greys")[2:5]) %>% 
  setNames(.,clindata$`Sample Name` %>% sort(decreasing=TRUE))