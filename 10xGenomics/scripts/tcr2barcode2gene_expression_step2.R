library(Seurat)
args= commandArgs(T)

RNAseq_path_10x = args[1]
# e.g. P0000-blackbird/2103-HC001-ZJW/G85E2L1_RNA-Seq/filtered_feature_bc_matrix
candidate_barcode = args[2] 
# e.g. "total.TCR.clonotype.in.raw_barcode.csv"
outputfile = args[3]

HC = Read10X(RNAseq_path_10x)
HC = CreateSeuratObject(counts = HC, min.cells = 3)
###################
# read TCR data
TCR_barcode_info = read.csv(candidate_barcode)
HC@meta.data = cbind(HC@meta.data, TCR_barcode_info[,-1])

#######################################################################################
# # normalization
# HC = NormalizeData(HC)
# # find variable features
# HC = FindVariableFeatures(HC)

# top10 = head(VariableFeatures(HC), 10)
# top10
# run pca
# HC = ScaleData(HC)
# HC = RunPCA(HC, features = VariableFeatures(object = HC))
# DimPlot(HC, reduction = "pca")  == PCAPlot(HC)

# HC <- FindNeighbors(HC, dims = 1:10)
# HC <- FindClusters(HC, resolution = 0.5)
# # Run non-linear dimensional reduction (UMAP/tSNE)
# HC <- RunUMAP(HC, dims = 1:10)
# HC <- RunTSNE(HC)
# DimPlot(HC, reduction = "umap")
# DimPlot(HC, reduction = "tsne")

# UMAPPlot(HC, group.by = "expansion")
# TSNEPlot(HC, group.by = "expansion")
###########################################################################################
cell.use = rownames(HC@meta.data[which(HC@meta.data$expansion =="Known"),])
sub_HC = subset(HC, cells =cell.use)
######### filtering MT
# sub_HC[["percent.mt"]] = PercentageFeatureSet(sub_HC, pattern = "^MT-")
# sub_HC = subset(sub_HC, subset = percent.mt<5)
# dim(sub_HC@meta.data)

# convert a sparse matrix to a data frame like this:
gene_expression_matrix = as.data.frame(sub_HC@assays$RNA@counts)
# filter gene by 0 line.
# gene_expression_matrix_filtered = gene_expression_matrix[which(rowSums(gene_expression_matrix)>0),]
gene_expression_matrix_filtered = gene_expression_matrix
# rbind barcode info to matrix. 
info.use = TCR_barcode_info[which(TCR_barcode_info$expansion =="Known"),]
info.use.t = t(info.use)
colnames(info.use.t) =  info.use.t[1,]

rbind_info = rbind(info.use.t, gene_expression_matrix_filtered, stringsAsFactors = FALSE)

# # useless
# gene_expression_matrix_filtered_combind_info = cbind(info.use[-4], gene_expression_matrix_filtered)
# gene_expression_matrix_filtered_combind_info[1:10,1:10]
library(tibble)
write.table(rbind_info[-c(1,4),] %>% rownames_to_column("barcode"), file = outputfile, sep = "\t", row.names =FALSE, quote = FALSE)
