## ---- add_clinical_features
seurat_obj@meta.data$sex <- clindata$Sex[match(seurat_obj@meta.data$orig.ident, clindata$`Sample Name`)]
seurat_obj@meta.data$age <- clindata$Age[match(seurat_obj@meta.data$orig.ident, clindata$`Sample Name`)]
seurat_obj@meta.data$Genotype <- clindata$Genotype[match(seurat_obj@meta.data$orig.ident, clindata$`Sample Name`)]
seurat_obj@meta.data$cell_id <- rownames(seurat_obj@meta.data)
