library(Seurat)
library(dplyr)
args = commandArgs(T)

cellranger_folder = args[1]
data = Read10X(data.dir = paste(cellranger_folder,'outs/filtered_feature_bc_matrix',sep="/"))

########## with CITEseq 
seurat_obj = CreateSeuratObject(data$`Gene Expression`)
seurat_obj = NormalizeData(seurat_obj) %>% FindVariableFeatures() %>% ScaleData() %>% 
  RunPCA() %>% RunUMAP(dims=1:50) %>% FindNeighbors(dims=1:50) %>% 
    FindClusters(resolution=0.2)

write.csv(Cells(seurat_obj), file = "cellID_obs.csv", row.names = FALSE)
write.csv(Embeddings(seurat_obj, reduction = "umap"), file = "cell_embeddings.csv")
write.csv(seurat_obj@meta.data[,'seurat_clusters',drop=F], file = "cell_clusters.csv")

