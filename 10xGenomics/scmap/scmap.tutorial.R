library(scmap)
browseVignettes("scmap")
# BiocManager::install("scmap")

# 
# corrplot function in R provides a nice heatmap-style visualization of this.
# 
# 
# CorOb.cor.exp <- as.data.frame(cor(CorOb.av.exp))
# CorOb.cor.exp$x <- rownames(CorOb.cor.exp)
# neupeptide.cor.df <- tidyr::gather(data = CorOb.cor.exp, y, correlation,CorOb.cor.exp$x)
# a <- ggplot(neupeptide.cor.df, aes(x, y, fill = correlation)) + geom_raster()+scale_fill_gradientn(colors = c("blue", "white", "red"))+ 
#   theme(axis.text.x = element_text(angle = 45, vjust = 1, 
#                                    size = 12, hjust = 1))+
#   coord_fixed()
# 
# Two scRNA-seq tools, namely SingleR and RCA (Reference Component Analysis, Li et. al, 2017) compute correlations of bulk transcriptomic 
# data to single cells; for this they use Spearman correlation and Pearson correlation respectively (The two tools use the same bulk 
#                                                                                                    transcriptome reference). 
# I am using these two tools for different purposes so I did not really compare their results, however, 
# single cell clusters that I compute using the reference components (calcuated with RCA) is in good agreement with single cells 
# labels computed with RCA (adjusted rand index of 0.81).
# 
# https://bioinformatics.stackexchange.com/questions/10582/which-correlation-method-to-compute-the-correlation-score-between-different-clus
# 

setwd("/cygene2/work/P0000-Blackbird/2103-GB001/GB001003/dsy_GB13_Tumor.Digest_hCD45_VS_GB13_TU_CL1_analysis")
list.files()
folder1 = "../GB13_TU_CL1_G0172E5L1_scRNAseq/outs/filtered_feature_bc_matrix"
folder2 = "../GB13_TumorDigest-hCD45_G0172E4L1_scRNAseq/outs/filtered_feature_bc_matrix"

library(Seurat)
library(SingleCellExperiment)
library(dplyr)
# read 10x data
data1 <- Seurat::Read10X(data.dir = folder1)
data2 <- Seurat::Read10X(data.dir = folder2)

seurat_obj1 = CreateSeuratObject(data1)
seurat_obj2 = CreateSeuratObject(data2)

# since we don't have cell type information, so just using cluster for annotation.
standard_seurat_pipeline = function(seurat_obj){
  seurat_obj= NormalizeData(seurat_obj) %>% FindVariableFeatures() %>% 
    ScaleData() %>% RunPCA() %>% FindNeighbors() %>% 
    FindClusters(resolution = 0.2) %>% RunUMAP(dims=1:50)
  return(seurat_obj)
}
seurat_obj1 = standard_seurat_pipeline(seurat_obj1)
seurat_obj2 = standard_seurat_pipeline(seurat_obj2)

DimPlot(seurat_obj1) |DimPlot(seurat_obj2)

# convert seurat object to SingleCellExperiment object.
sce_1 <- Seurat::as.SingleCellExperiment(seurat_obj1)
sce_2 <- Seurat::as.SingleCellExperiment(seurat_obj2)

rowData(sce_1)$feature_symbol <- rownames(sce_1)
rowData(sce_2)$feature_symbol <- rownames(sce_2)
sce_1 <- sce_1[!duplicated(rownames(sce_1)), ]
sce_2 <- sce_2[!duplicated(rownames(sce_2)), ]

sce_1 <- selectFeatures(sce_1, suppress_plot = FALSE)
sce_2 <- selectFeatures(sce_2, suppress_plot = FALSE)


# scmap-cluster default: cell_type1
sce_2 <- indexCluster(sce_2, cluster_col = "seurat_clusters" )
head(metadata(sce_2)$scmap_cluster_index)
heatmap(as.matrix(metadata(sce_2)$scmap_cluster_index))

############## projection
scmapCluster_results <- scmapCluster(
  projection = sce_1,
  index_list = list(yan = metadata(sce_2)$scmap_cluster_index)
)
### results
head(scmapCluster_results$scmap_cluster_labs)
head(scmapCluster_results$scmap_cluster_siml)
head(scmapCluster_results$combined_labs)

plot(
  getSankey(
    colData(sce_1)$seurat_clusters,
    scmapCluster_results$scmap_cluster_labs[,'yan'],
    plot_height = 400
  )
)
