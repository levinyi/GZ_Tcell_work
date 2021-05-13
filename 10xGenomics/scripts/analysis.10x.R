library(Seurat)
library(tidyverse)
args = commandArgs(T)

## setup for windows
# setwd("C:\\Users\\dushiyi\\Nutstore\\.nutstore_c2hpeWlAcm9vdHBhdGhneC5jb20=\\DuShiYi\\P0000-blackbird\\2103-CR001\\CR001004\\CR001004_CD4_CD8_scRNA_TCRseq")
# list.dirs()
# scRNAseq_path = "CR001004_G179E4L1_scRNAseq"
# tcr_path = "CR001004_G179E3L1_TCRseq_IMGT"
# # something for linux
scRNAseq_path = args[1]
tcr_path = args[2]

# read scRNAseq folder:
data = Seurat::Read10X(data.dir = paste(scRNAseq_path,'filtered_feature_bc_matrix',sep = "/"))
sc_seurat_obj <- Seurat::CreateSeuratObject(data, project = "scRNAseq_object")
# sc_seurat_obj # 33538 features across 5533 samples within 1 assay 

# read TCRseq file
add_clonotype <- function(tcr_path, seurat_obj){
  tcr_data = read.csv(paste(tcr_path, "filtered_contig_annotations.csv", sep = "/"))
  tcr_data <- tcr_data[!duplicated(tcr_data$barcode), ] %>% 
    dplyr::select(barcode, raw_clonotype_id) %>% 
    dplyr::rename(clonotype_id = raw_clonotype_id)
  # Clonotype-centric info.
  clono <- read.csv(paste(tcr_path,"clonotypes.csv", sep="/"))
  head(clono[, c("clonotype_id", "cdr3s_aa")])
  # Slap the AA sequences onto our original table by clonotype_id.
  tcr_data <- merge(tcr_data, clono[, c("clonotype_id", "cdr3s_aa")])
  # Reorder so barcodes are first column and set them as rownames.
  tcr_data <- tcr_data %>% tibble::column_to_rownames("barcode")

  # Add to the Seurat object's metadata.
  clono_seurat <- Seurat::AddMetaData(object=sc_seurat_obj, metadata=tcr_data)
  return(clono_seurat)
}
sc_seurat_obj = add_clonotype(tcr_path, sc_seurat_obj)


################
sc_seurat_obj[["percent.mt"]] <- PercentageFeatureSet(sc_seurat_obj, pattern = "^MT-")
p1 = VlnPlot(sc_seurat_obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
ggsave(filename="P1.VlnPlot.feature.png", plot=p1, path=".")
# FeatureScatter is typically used to visualize feature-feature relationships, but can be used
# for anything calculated by the object, i.e. columns in object metadata, PC scores etc.

plot1 <- FeatureScatter(sc_seurat_obj, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(sc_seurat_obj, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
p2 = CombinePlots(plots = list(plot1, plot2))
ggsave(filename="P2.combinePlots.FeatureScatter.png", plot=p2, path=".")

sc_seurat_obj <- subset(sc_seurat_obj, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)
# sc_seurat_obj # 33538 features across 3549 samples within 1 assay 
# standard seurat workflow
sc_seurat_obj = NormalizeData(object = sc_seurat_obj)
sc_seurat_obj = FindVariableFeatures(object = sc_seurat_obj)
sc_seurat_obj = ScaleData(object = sc_seurat_obj)
sc_seurat_obj = RunPCA(object = sc_seurat_obj, npcs = 20)
sc_seurat_obj = FindNeighbors(object = sc_seurat_obj)
sc_seurat_obj = FindClusters(object = sc_seurat_obj, resolution = 0.5)
# sc_seurat_obj = RunTSNE(object = sc_seurat_obj)
sc_seurat_obj = RunUMAP(object = sc_seurat_obj, dims = 1:20)
p3 <- DimPlot(sc_seurat_obj, reduction = "umap")
ggsave(filename = "P3.raw.umap.png",plot = p3, path = "." )

# find doublets
## doubletfinder
library(DoubletFinder)
find_doublets <- function(sc_seurat_obj){
  sweep.res.list <- DoubletFinder::paramSweep_v3(sc_seurat_obj, PCs = 1:10, sct = FALSE)
  sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)
  bcmvn <- find.pK(sweep.stats)
  mpK <- as.numeric(as.vector(bcmvn$pK[which.max(bcmvn$BCmetric)]))
  ## ----------------- DoubletFinder:Homotypic Doublet Proportion Estimate ------------------------
  annotations <- sc_seurat_obj@meta.data$seurat_clusters
  homotypic.prop <- modelHomotypic(annotations) 
  nExp_poi <- round(0.075*nrow(sc_seurat_obj@meta.data)) ## Assuming 7.5% doublet formation rate
  nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
  
  ## ------------DoubletFinder:Run DoubletFinder with varying classification stringencies ---------
  sc_seurat_obj <- doubletFinder_v3(sc_seurat_obj, PCs = 1:10, pN = 0.25, pK = mpK, 
                                    nExp = nExp_poi, reuse.pANN = FALSE, sct = TRUE)
  sc_seurat_obj <- doubletFinder_v3(sc_seurat_obj, PCs = 1:10, pN = 0.25, pK = mpK, 
                                    nExp = nExp_poi.adj, reuse.pANN = paste("pANN_0.25",mpK,nExp_poi,sep = "_"), sct = TRUE)
  
  ## ---------------------Doublets by DoubletFinder --------------------------------------------
  sc_seurat_obj$DF = sc_seurat_obj@meta.data[,match(paste("DF.classifications_0.25",mpK,nExp_poi,sep = "_"),colnames(sc_seurat_obj@meta.data))]
  sc_seurat_obj$DF[which(sc_seurat_obj$DF == "Doublet" && sc_seurat_obj@meta.data[,match(paste("DF.classifications_0.25",mpK,nExp_poi.adj,sep = "_"),colnames(sc_seurat_obj@meta.data))])] <- "Doublet_lo"
  sc_seurat_obj$DF[which(sc_seurat_obj$DF == "Doublet")] <- "Doublet_hi"
  sc_seurat_obj$DF
  write.table(sc_seurat_obj$DF, file = "P4.doublets.table.csv", sep = ",", row.names = FALSE, quote = FALSE)
  p4 = DimPlot(sc_seurat_obj, reduction = "umap", group.by = "DF")
  ggsave(filename = "P4.doublets.png",plot = p4, path = "./" )
  return(sc_seurat_obj)
}
sc_seurat_obj = find_doublets(sc_seurat_obj)
sc_seurat_obj = subset(sc_seurat_obj,subset=DF=='Singlet')
# sc_seurat_obj # 33538 features across 3283 samples within 1 assay 
# singR annotate cluster
library(SingleR)
library(celldex)
cluster_annotation <- function(sc_seurat_obj){
  sc_seurat_obj_SingleR = GetAssayData(sc_seurat_obj, slot="data")
  sc_seurat_obj_clusters = sc_seurat_obj@meta.data$seurat_clusters
  
  hpca = celldex::HumanPrimaryCellAtlasData()
  blueprint_encode = celldex::BlueprintEncodeData()
  dice = celldex::DatabaseImmuneCellExpressionData()
  mona = celldex::MonacoImmuneData()
  novershter = celldex::NovershternHematopoieticData()
  
  pred.hpca = SingleR(test = sc_seurat_obj_SingleR, ref = hpca, labels = hpca$label.main,clusters = sc_seurat_obj_clusters )
  pred.blue = SingleR(test = sc_seurat_obj_SingleR, ref = blueprint_encode,labels = blueprint_encode$label.main,clusters = sc_seurat_obj_clusters )
  pred.dice = SingleR(test = sc_seurat_obj_SingleR, ref =dice, labels = dice$label.main,clusters = sc_seurat_obj_clusters)
  pred.mona = SingleR(test = sc_seurat_obj_SingleR, ref = mona, labels = mona$label.main,clusters = sc_seurat_obj_clusters)
  pred.novershter = SingleR(test = sc_seurat_obj_SingleR, ref = novershter, labels = novershter$label.main,clusters = sc_seurat_obj_clusters)
  
  sc_seurat_obj_cellType = data.frame(ClusterID=levels(sc_seurat_obj@meta.data$seurat_clusters),
                                      hpca = pred.hpca$labels,
                                      blue  = pred.blue$labels,
                                      Dice = pred.dice$labels,
                                      mona = pred.mona$labels,
                                      novershter = pred.novershter$labels
  )
  write.table(sc_seurat_obj_cellType, file = "P5.SingleR_cell_type.csv",sep = ",", row.names = FALSE, quote=FALSE)
  sc_seurat_obj@meta.data$singleR = sc_seurat_obj_cellType[match(sc_seurat_obj_clusters, sc_seurat_obj_cellType$ClusterID),'hpca']
  p5_0 = DimPlot(sc_seurat_obj, reduction = "umap", label = T, group.by = 'singleR',)+ggtitle("HumanPrimaryCellAtlasData")
  sc_seurat_obj@meta.data$singleR = sc_seurat_obj_cellType[match(sc_seurat_obj_clusters, sc_seurat_obj_cellType$ClusterID),'blue']
  p5_1 = DimPlot(sc_seurat_obj, reduction = "umap", label = T, group.by = 'singleR',)+ggtitle("BlueprintEncodeData")
  sc_seurat_obj@meta.data$singleR = sc_seurat_obj_cellType[match(sc_seurat_obj_clusters, sc_seurat_obj_cellType$ClusterID),'Dice']
  p5_2 = DimPlot(sc_seurat_obj, reduction = "umap", label = T, group.by = 'singleR',)+ggtitle("DatabaseImmuneCellExpressionData")
  sc_seurat_obj@meta.data$singleR = sc_seurat_obj_cellType[match(sc_seurat_obj_clusters, sc_seurat_obj_cellType$ClusterID),'mona']
  p5_3 = DimPlot(sc_seurat_obj, reduction = "umap", label = T, group.by = 'singleR') + ggtitle("MonacoImmuneData")
  sc_seurat_obj@meta.data$singleR = sc_seurat_obj_cellType[match(sc_seurat_obj_clusters, sc_seurat_obj_cellType$ClusterID),'novershter']
  p5_4 = DimPlot(sc_seurat_obj, reduction = "umap", label = T, group.by = 'singleR',)+ ggtitle("NovershternHematopoieticData")
  ggsave(filename="P5_0.cluster.annotation.with.hpca.png",plot = p5_0, path = ".")
  ggsave(filename="P5_1.cluster.annotation.with.blue.png",plot = p5_1, path = ".")
  ggsave(filename="P5_2.cluster.annotation.with.dice.png",plot = p5_2, path = ".")
  ggsave(filename="P5_3.cluster.annotation.with.mona.png",plot = p5_3, path = ".")
  ggsave(filename="P5_4.cluster.annotation.with.nove.png",plot = p5_4, path = ".")
  # library(cowplot)
  # p5 = cowplot::plot_grid(p5_0,p5_1,p5_2,p5_3,p5_4, ncol = 2)
  # ggsave(filename = "cluster.annotation.with.5.databases.png", plot = p5, path = "./" )
  return(sc_seurat_obj)
}
sc_seurat_obj = cluster_annotation(sc_seurat_obj)
save.image("sc_seurat_obj")

