## ---- calculate_percent_mt
seurat_obj <- PercentageFeatureSet(seurat_obj,
                                   pattern = "^MT-",
                                   col.name = "percent_mt")

## ---- find_variable_features
seurat_obj <- FindVariableFeatures(seurat_obj)

## ---- calculate_doublet
set.seed(10)
doubletDensity <- doubletCells(seurat_obj@assays$RNA@counts,
                               subset.row=seurat_obj@assays$RNA@meta.features$vst.variable,
                               d=max(seuratAnalysisCfg$PCA_DIMS))
seurat_obj@meta.data$doublet_score <- log10(doubletDensity+1)

## ---- filter_cell
i <- seurat_obj@meta.data$nFeature_RNA >= seuratAnalysisCfg$MIN_FEATURES
j <- seurat_obj@meta.data$doublet_score < seuratAnalysisCfg$MAX_DOUBLET_SCORE
k <- seurat_obj@meta.data$percent_mt < seuratAnalysisCfg$PERC_MITO
seurat_obj <- seurat_obj[, i&j&k]
## ---- filter_genes_low
#remove genes detected in only a 20 cells
seurat_obj <- seurat_obj[rowSums(seurat_obj@assays$RNA@counts>0) >= seuratAnalysisCfg$MIN_CELLS ,]


