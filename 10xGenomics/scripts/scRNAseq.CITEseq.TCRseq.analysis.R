# install.packages("Seurat")
# install.packages("tidyr")
# install.packages("tidyverse")
library(Seurat)
library(sctransform)
library(tidyverse, quietly=T)

args = commandArgs(T)

# for 3 ADT:
# adt_CD4-TotalSeqC, adt_CD8-TotalSeqC, adt_IgG1-TotalSeqC

## setup for windows
# setwd("C:\\Users\\CEID\\Nutstore\\.nutstore_c2hpeWlAcm9vdHBhdGhneC5jb20=\\DuShiYi\\P0000-blackbird\\2103-CR001\\CR001003\\CR001003_PD1_CD4_CD8\\G082E5L1_scRNAseq_G082E6L1_CITEseq")
# list.dirs()
scRNAseq_path = "./"
tcr_path = "../../G471E1L1_TCR_IMGT"

## setting for linux
scRNAseq_path = args[1]   # scRNAseq data must contain ADT information.
tcr_path = args[2]


# create a output directory in current path.
dir_list <- c("scRNAseq_CITEseq_TCRseq_analysis","scRNAseq_CITEseq_TCRseq_analysis/1QC",
              "scRNAseq_CITEseq_TCRseq_analysis/2Cluster","scRNAseq_CITEseq_TCRseq_analysis/3CellMarker",
              "scRNAseq_CITEseq_TCRseq_analysis/4Annotation","scRNAseq_CITEseq_TCRseq_analysis/5TrajectoryAnalysis",
              "scRNAseq_CITEseq_TCRseq_analysis/6RNAvelocity","scRNAseq_CITEseq_TCRseq_analysis/7CellPhone",
              "scRNAseq_CITEseq_TCRseq_analysis/4Annotation/SingleR","scRNAseq_CITEseq_TCRseq_analysis/4Annotation/Garnett",
              "scRNAseq_CITEseq_TCRseq_analysis/4Annotation/scCATCH","scRNAseq_CITEseq_TCRseq_analysis/",
              "scRNAseq_CITEseq_TCRseq_analysis/")
for (each in dir_list){
  if(!dir.exists(each)){
    dir.create(each)
  }
}

# dir.create("scRNAseq_CITEseq_TCRseq_analysis")
output_dir = "scRNAseq_CITEseq_TCRseq_analysis"

# read scRNAseq folder:
data = Seurat::Read10X(data.dir = paste(scRNAseq_path, '/filtered_feature_bc_matrix', sep = "/"))
# 10X data contains more than one type and is being returned as a list containing matrices of each type.

output_ADT_data = data$`Antibody Capture` %>% as.matrix() %>% t() %>% as.data.frame() %>% rownames_to_column("barcode")
write.table(output_ADT_data, file=paste(output_dir,'ADT.t.csv',sep="/"), sep = ",",row.names = FALSE)

output_RNA_data = data$`Gene Expression` %>% as.matrix() %>% t() %>%  as.data.frame() %>%   rownames_to_column("barcode") %>% 
  select(c("barcode","CD4","CD8A","CD8B","CD3E"))
write.table(output_RNA_data, file=paste(output_dir,'RNA.t.csv',sep="/"), sep = ",", row.names = FALSE)

sc_seurat_obj <- CreateSeuratObject(data$`Gene Expression`, project = "scRNAseq_object") 

# Assays(sc_seurat_obj)
# [1] "RNA" "SCT"
# > sc_seurat_obj
# An object of class Seurat 
# 52716 features across 5271 samples within 2 assays 
# Active assay: SCT (16115 features, 3000 variable features)
# 1 other assay present: RNA
# 2 dimensional reductions calculated: pca, umap

sc_seurat_obj[['ADT']] <- CreateAssayObject(counts = data$`Antibody Capture`)

# Validate that the object now contains multiple assays
Assays(sc_seurat_obj)
## [1] "RNA" "ADT" "ADT"
rownames(sc_seurat_obj[["ADT"]])
## [1] "CD8-TotalSeqC" "CD4-TotalSeqC" "IgG1-TotalSeqC"


# read TCRseq file
add_clonotype <- function(tcr_path, seurat_obj){
  tcr_data = read.csv(paste(tcr_path, "outs/filtered_contig_annotations.csv", sep = "/"))
  tcr_data <- tcr_data[!duplicated(tcr_data$barcode), ] %>% 
    dplyr::select(barcode, raw_clonotype_id) %>% 
    dplyr::rename(clonotype_id = raw_clonotype_id)
  # Clonotype-centric info.
  clono <- read.csv(paste(tcr_path,"outs/clonotypes.csv", sep="/"))
  # head(clono[, c("clonotype_id", "cdr3s_aa")])
  # Slap the AA sequences onto our original table by clonotype_id.
  tcr_data <- merge(tcr_data, clono[, c("clonotype_id", "cdr3s_aa")])
  # Reorder so barcodes are first column and set them as rownames.
  tcr_data <- tcr_data %>% tibble::column_to_rownames("barcode")

  # Add to the Seurat object's metadata.
  clono_seurat <- Seurat::AddMetaData(object=seurat_obj, metadata=tcr_data)
  return(clono_seurat)
}
sc_seurat_obj = add_clonotype(tcr_path, sc_seurat_obj)


################
sc_seurat_obj = PercentageFeatureSet(sc_seurat_obj, pattern = "^MT-", col.name = "percent.mt")
# log scale, flip=FALSE is default.
p1 = VlnPlot(sc_seurat_obj, assay = "RNA", features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, log=TRUE, flip=FALSE)
p1
ggsave(filename=paste(output_dir,"1QC/P1.VlnPlot.feature.png",sep="/"), plot=p1, width=1344,height=960, units="px", path=".")

# FeatureScatter is typically used to visualize feature-feature relationships, but can be used
# for anything calculated by the object, i.e. columns in object metadata, PC scores etc.

plot1 <- FeatureScatter(sc_seurat_obj, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(sc_seurat_obj, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
p2 = CombinePlots(plots = list(plot1, plot2))
p2
ggsave(filename=paste(output_dir,"1QC/P2.combinePlots.FeatureScatter.png",sep="/"), plot=p2, width=1344,height=960, units = "px", path=".")

# sc_seurat_obj <- subset(sc_seurat_obj, subset = nFeature_RNA > 20 & nFeature_RNA < 3000 & percent.mt < 25)

sc_seurat_obj = SCTransform(sc_seurat_obj, vars.to.regress = "percent.mt")
sc_seurat_obj = RunPCA(sc_seurat_obj, npcs = 50)
######
VizDimLoadings(sc_seurat_obj, dims = 1:2, reduction = "pca")
DimPlot(sc_seurat_obj, reduction = "pca")
DimHeatmap(sc_seurat_obj, dims = 1, cells = 500, balanced = TRUE)
DimHeatmap(sc_seurat_obj, dims = 1:15, cells = 500, balanced = TRUE)

# Assays(sc_seurat_obj) <- "RNA"
# DefaultAssay(sc_seurat_obj) <- "RNA"
# # Determin the "dimensionality" of the dataset.
# sc_seurat_obj <- JackStraw(sc_seurat_obj, num.replicate = 100)
# sc_seurat_obj <- ScoreJackStraw(sc_seurat_obj, dims = 1:20)
ElbowPlot(sc_seurat_obj,ndims = 50)

sc_seurat_obj = FindNeighbors(sc_seurat_obj)
sc_seurat_obj = RunUMAP(sc_seurat_obj,dims = 1:20)
sc_seurat_obj = RunTSNE(sc_seurat_obj)
sc_seurat_obj = FindClusters(sc_seurat_obj, resolution = 0.8)


# sc_seurat_obj # 33538 features across 3549 samples within 1 assay 
# standard seurat workflow
#sc_seurat_obj = RunPCA(object = sc_seurat_obj, npcs = 20)
#sc_seurat_obj = FindNeighbors(object = sc_seurat_obj)
#sc_seurat_obj = FindClusters(object = sc_seurat_obj, resolution = 0.5)
#sc_seurat_obj = RunTSNE(object = sc_seurat_obj)
#sc_seurat_obj = RunUMAP(object = sc_seurat_obj, dims = 1:20)

p3 <- DimPlot(sc_seurat_obj, reduction = "umap")
p3_2 <- DimPlot(sc_seurat_obj, reduction = "tsne")
p3
p3_2
ggsave(filename = paste(output_dir,"/2Cluster/P3.raw.umap.png",sep="/"),plot = p3, width=1344,height=960,units="px",path = "." )
ggsave(filename = paste(output_dir,"/2Cluster/P3.raw.tsne.png",sep="/"),plot = p3_2, width=1344,height=960,units="px",path = "." )


########################################
########################################
# for RNA velocity
# save cells to RNA velocity analysis
write.csv(Seurat::Cells(sc_seurat_obj), file = paste(output_dir,"6RNAvelocity/cellID_obs.csv", sep = "/"), row.names = FALSE)
write.csv(Seurat::Embeddings(sc_seurat_obj, reduction = "umap"), file = paste(output_dir,"6RNAvelocity/cell_embeddings.csv", sep = "/"))
write.csv(sc_seurat_obj@meta.data$seurat_clusters, file = paste(output_dir,"6RNAvelocity/clusters.csv", sep = "/"))
########################################
########################################
# for PieParty: https://github.com/harbourlab/PieParty
write.csv(Seurat::GetAssayData(sc_seurat_obj, slot = "counts"),
	  file = paste(output_dir,"6RNAvelocity/expression.counts.for.PieParty.csv", sep = "/"), 
	  row.names = FALSE)
#########################################
#########################################
# find doublets
## doubletfinder
## install.packages('remotes')
## remotes::install_github('chris-mcginnis-ucsf/DoubletFinder')
library(DoubletFinder)

# Do not apply DoubletFinder to aggregated scRNA-seq data representing multiple distinct samples
## -------------- pK Identification (no ground-truth) --------------------------------------
sweep.res.list <- DoubletFinder::paramSweep_v3(sc_seurat_obj, PCs = 1:20, sct = FALSE)
# sweep.res.list <- DoubletFinder::paramSweep_v3(sc_seurat_obj, PCs = 1:15, sct = T) # need test
sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)
bcmvn <- find.pK(sweep.stats)
mpK <- as.numeric(as.vector(bcmvn$pK[which.max(bcmvn$BCmetric)]))

## ----------------- DoubletFinder:Homotypic Doublet Proportion Estimate ------------------------
annotations <- sc_seurat_obj@meta.data$seurat_clusters
homotypic.prop <- modelHomotypic(annotations) 
nExp_poi <- round(0.075*nrow(sc_seurat_obj@meta.data)) ## Assuming 7.5% doublet formation rate
# DoubletRate = 0.039  # 5000 Cells correspond to doublets rate yes 3.9%
# nExp_poi <- round(DoubletRate*ncol(sc_seurat_obj))
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

## ------------DoubletFinder:Run DoubletFinder with varying classification stringencies ---------
sc_seurat_obj <- doubletFinder_v3(sc_seurat_obj, PCs = 1:10, pN = 0.25, pK = mpK, 
                                  nExp = nExp_poi, reuse.pANN = FALSE, sct = TRUE)
sc_seurat_obj <- doubletFinder_v3(sc_seurat_obj, PCs = 1:10, pN = 0.25, pK = mpK, 
                                  nExp = nExp_poi.adj, reuse.pANN = paste("pANN_0.25",mpK,nExp_poi,sep = "_"), sct = TRUE)

## ---------------------Doublets by DoubletFinder --------------------------------------------
sc_seurat_obj$DF = sc_seurat_obj@meta.data[,match(paste("DF.classifications_0.25",mpK,nExp_poi,sep = "_"),colnames(sc_seurat_obj@meta.data))]
doublets_rate = round(length(which(sc_seurat_obj$DF=="Doublet"))/ nrow(sc_seurat_obj@meta.data)*100,2)
print(paste("doublets_rate : ", doublets_rate,"%", sep=""))
sc_seurat_obj$DF[which((sc_seurat_obj$DF == "Doublet") & (sc_seurat_obj@meta.data[,match(paste("DF.classifications_0.25",mpK,nExp_poi.adj,sep = "_"),colnames(sc_seurat_obj@meta.data))] == "Singlet"))] <- "Doublet_lo"
sc_seurat_obj$DF[which(sc_seurat_obj$DF == "Doublet")] <- "Doublet_hi"

p4 = DimPlot(sc_seurat_obj, reduction = "umap", group.by = "DF", label.size = 5, label = TRUE, pt.size = 1) + 
  ggtitle(paste("Doublets rate : ",doublets_rate,"%", sep=""))
ggsave(filename = paste(output_dir, "1QC/P4.doublets.png",sep="/"), plot = p4, width=9, height=7, path = "./" )

########################## finish Doubletsfinder.
##########################################################
#########################################################

## select Singlet cells 
sc_seurat_obj = subset(sc_seurat_obj, subset=DF=='Singlet')
## sc_seurat_obj # 33538 features across 3283 samples within 1 assay 
#########################################################
#########################################################
############################# singR annotate cluster
## install.packages("BiocManager")
## BiocManager::install("SingleR")
## or 
## install.packages("devtools")
## devtools::install_github('dviraran/SingleR')

library(SingleR)
# BiocManager::install("celldex")
library(celldex)

sc_seurat_obj_SingleR = GetAssayData(sc_seurat_obj, slot="data")
sc_seurat_obj_clusters = sc_seurat_obj@meta.data$seurat_clusters

hpca = celldex::HumanPrimaryCellAtlasData()
blueprint_encode = celldex::BlueprintEncodeData()
dice = celldex::DatabaseImmuneCellExpressionData()
mona = celldex::MonacoImmuneData()
novershter = celldex::NovershternHematopoieticData()

# choose label.main or label.fine
pred.hpca = SingleR(test = sc_seurat_obj_SingleR, ref = hpca, labels = hpca$label.fine, clusters = sc_seurat_obj_clusters )
pred.blue = SingleR(test = sc_seurat_obj_SingleR, ref = blueprint_encode,labels = blueprint_encode$label.fine,clusters = sc_seurat_obj_clusters )
pred.dice = SingleR(test = sc_seurat_obj_SingleR, ref =dice, labels = dice$label.fine,clusters = sc_seurat_obj_clusters)
pred.mona = SingleR(test = sc_seurat_obj_SingleR, ref = mona, labels = mona$label.fine,clusters = sc_seurat_obj_clusters)
pred.novershter = SingleR(test = sc_seurat_obj_SingleR, ref = novershter, labels = novershter$label.fine,clusters = sc_seurat_obj_clusters)

sc_seurat_obj_cellType = data.frame(ClusterID=levels(sc_seurat_obj@meta.data$seurat_clusters),
                                    hpca = pred.hpca$labels,
                                    blue  = pred.blue$labels,
                                    Dice = pred.dice$labels,
                                    mona = pred.mona$labels,
                                    novershter = pred.novershter$labels
)
write.table(sc_seurat_obj_cellType, file = paste(output_dir,"4Annotation/SingleR/P5.SingleR_cell_type.csv",sep="/"),sep = ",", row.names = FALSE, quote=FALSE)
sc_seurat_obj@meta.data$singleR.hpca = sc_seurat_obj_cellType[match(sc_seurat_obj_clusters, sc_seurat_obj_cellType$ClusterID),'hpca']
p5_0 = DimPlot(sc_seurat_obj, reduction = "umap", label = T, group.by = 'singleR.hpca',)+ggtitle("HumanPrimaryCellAtlasData")
sc_seurat_obj@meta.data$singleR.blue = sc_seurat_obj_cellType[match(sc_seurat_obj_clusters, sc_seurat_obj_cellType$ClusterID),'blue']
p5_1 = DimPlot(sc_seurat_obj, reduction = "umap", label = T, group.by = 'singleR.blue',)+ggtitle("BlueprintEncodeData")
sc_seurat_obj@meta.data$singleR.dice = sc_seurat_obj_cellType[match(sc_seurat_obj_clusters, sc_seurat_obj_cellType$ClusterID),'Dice']
p5_2 = DimPlot(sc_seurat_obj, reduction = "umap", label = T, group.by = 'singleR.dice',)+ggtitle("DatabaseImmuneCellExpressionData")
sc_seurat_obj@meta.data$singleR.mona = sc_seurat_obj_cellType[match(sc_seurat_obj_clusters, sc_seurat_obj_cellType$ClusterID),'mona']
p5_3 = DimPlot(sc_seurat_obj, reduction = "umap", label = T, group.by = 'singleR.mona') + ggtitle("MonacoImmuneData")
sc_seurat_obj@meta.data$singleR.nove = sc_seurat_obj_cellType[match(sc_seurat_obj_clusters, sc_seurat_obj_cellType$ClusterID),'novershter']
p5_4 = DimPlot(sc_seurat_obj, reduction = "umap", label = T, group.by = 'singleR.nove',)+ ggtitle("NovershternHematopoieticData")
ggsave(filename=paste(output_dir,"4Annotation/SingleR/P5_0.cluster.annotation.with.hpca.png",sep="/"),width=1344,height=960,units="px",plot = p5_0, path = "./")
ggsave(filename=paste(output_dir,"4Annotation/SingleR/P5_1.cluster.annotation.with.blue.png",sep="/"),width=1344,height=960,units="px",plot = p5_1, path = "./")
ggsave(filename=paste(output_dir,"4Annotation/SingleR/P5_2.cluster.annotation.with.dice.png",sep="/"),width=1344,height=960,units="px",plot = p5_2, path = "./")
ggsave(filename=paste(output_dir,"4Annotation/SingleR/P5_3.cluster.annotation.with.mona.png",sep="/"),width=1344,height=960,units="px",plot = p5_3, path = "./")
ggsave(filename=paste(output_dir,"4Annotation/SingleR/P5_4.cluster.annotation.with.nove.png",sep="/"),width=1344,height=960,units="px",plot = p5_4, path = "./")
library(cowplot)
p5 = cowplot::plot_grid(p5_0,p5_1,p5_2,p5_3,p5_4, ncol = 2)
# p5 = cowplot::plot_grid(p5_0,p5_1,p5_3,p5_4, ncol = 2)
ggsave(filename=paste(output_dir,"4Annotation/SingleR/P5_cluster.annotation.with.5.databases.png",sep="/"), plot = p5, width=1344,height=960,units="px", path = "./" )

########################################
########################################
######## scCATCH
library(scCATCH)
scCATCH_obj <- scCATCH::createscCATCH(sc_seurat_obj[["RNA"]]@data, cluster = as.character(sc_seurat_obj@meta.data$seurat_clusters))
clu_markers <- scCATCH::findmarkergene(scCATCH_obj, species = "Human", 
                                        cluster = "All", cancer = "Hepatocellular Cancer", 
                                        marker = cellmatch,
                                        tissue = c("Liver"), 
                                        use_method = "1", cell_min_pct = 0.25,
                                        logfc = 0.25, pvalue = 0.25,
                                        verbose = TRUE)
cellmatch.table = cellmatch
write.csv(cellmatch.table, file = paste(output_dir,  "4Annotation/scCATCH/CellMatch.database.csv", sep = "/"),row.names = FALSE, quote = FALSE)
clu_ann <- scCATCH::findcelltype(clu_markers)

scCATCH_ann_data = merge(clu_ann@meta, clu_ann@celltype,by = "cluster") %>% tibble::column_to_rownames("cell")
sc_seurat_obj = AddMetaData(sc_seurat_obj, scCATCH_ann_data)

p5_6 = DimPlot(sc_seurat_obj, reduction = "umap", label = T, group.by="cell_type", pt.size = 0.5) + ggtitle("scCATCH Peripheral blood")
p5_6
ggsave(filename=paste(output_dir, "4Annotation/scCATCH/P5_6.cluster.annotation.with.scCATCH.cellmatch.png", sep="/"), plot = p5_6, width=1344,height=960,units="px",path = "./" )

########################################
########################################
########  Garnett
library(monocle3)
library(garnett)
library(org.Hs.eg.db)
# Garnett 是基于monocle3， 所以它输入的数据格式是CellDataSet(CDS)
# create CDS object:
expression_matrix = GetAssayData(sc_seurat_obj, assay="RNA", slot = 'counts')
cell_metadata <- sc_seurat_obj@meta.data
gene_annotation <- data.frame(gene_short_name = rownames(expression_matrix))
rownames(gene_annotation) <- rownames(expression_matrix)

cds <- monocle3::new_cell_data_set(expression_matrix, cell_metadata = cell_metadata, gene_metadata = gene_annotation)
cds <- monocle3::preprocess_cds(cds, num_dim = 20)
# 接着做monocle3
p = plot_pc_variance_explained(cds)
p
ggsave(filename = paste(output_dir, "5TrajectoryAnalysis/P6.monocle3.pc_variance_explained.png", sep="/"), plot = p, width=1344,height=960,units="px", path = "./" )

cds <- monocle3::reduce_dimension(cds)
# pplot_cells(cds)
# head(rownames(cds))


cds <- cluster_cells(cds, reduction_method = "UMAP")
print("Learning graph, which can take a while depends on the sample")
cds <- learn_graph(cds, use_partition = T)

p = plot_cells(cds, genes=c("CD4","CD8A","CD8B"))
p
ggsave(filename = paste(output_dir, "5TrajectoryAnalysis/P6.monocle3.cells.CD4.CD8.png", sep="/"), plot = p, width=1344,height=960,units="px", path = "./" )

p = plot_cells(cds, color_cells_by = "seurat_clusters", label_cell_groups = F, )
p
ggsave(filename = paste(output_dir, "5TrajectoryAnalysis/P6.monocle3.color_cells_by_seurat_clusters.png", sep="/"), plot = p, width=1344,height=960,units="px", path = "./" )
### Construct and assign the made up partition

print("Plotting clusters")
library(ggrastr)
p = plot_cells(cds, x=1,y=2, 
           color_cells_by = "seurat_clusters",
           group_cells_by = "cluster", 
           # genes = c("CD4","CD8A","CD8B"),
           show_trajectory_graph = T,
           label_roots = F,
           graph_label_size = 3, 
           cell_size = .5,
           rasterize = T,
           alpha = .2,
           scale_to_range = 1,
           # label_principal_points = 1,
           label_groups_by_cluster = FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE,
)
p
ggsave(filename = paste(output_dir, "5TrajectoryAnalysis/P6.monocle3.cluster.png", sep="/"), plot = p, width=1344,height=960,units="px", path = "./" )


# #############
cds <- order_cells(cds, reduction_method = "UMAP")

print("Plotting color cell by pseudotime")
p = plot_cells(cds, x=1,y=2, 
           color_cells_by = 'pseudotime',
           group_cells_by = "cluster",
           # genes = c("CD4","CD8A","CD8B"),
           show_trajectory_graph = T,
           label_roots = F,
           graph_label_size = 3, 
           cell_size = .5,
           rasterize = T,
           alpha = .2,
           scale_to_range = 1,
           # label_principal_points = 1,
           label_groups_by_cluster = FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE,
           )
p
ggsave(filename = paste(output_dir, "5TrajectoryAnalysis/P6.monocle3.pseudotime.png", sep="/"), plot = p, width=1344,height=960,units="px", path = "./" )

# 单个基因的逆时序轨迹
AFD_genes = c("S100A9", "RPS27", "FCER1A")
AFD_lineage_cds = cds[AFD_genes,
                      clusters(cds) %in% c(11, 12, 5)]

p = plot_genes_in_pseudotime(AFD_lineage_cds,
                         color_cells_by="pseudotime",
                         min_expr=0.5)

ggsave(filename = paste(output_dir, "P6.monocle3.genes.in.pseudotime.png", sep="/"), plot = p, width=1344,height=960,units="px", path = "./" )

# garnett
##############################
# 2.2 marker 文件准备
# download.file(url="https://cole-trapnell-lab.github.io/garnett/marker_files/hsPBMC_markers.txt", destfile = "hsPBMC_markers.txt")
# 演示利用marker file训练分类器
# 2.3 marker 基因评估
# 对marker file中的marker基因评分
marker_file = "/cygene/work/00.test/pipeline/10xGenomics/scripts/RootPathPBMC_TCell_markers.txt"
marker_check <- check_markers(cds, marker_file,
                              db=org.Hs.eg.db,
                              cds_gene_id_type = "SYMBOL",
                              marker_file_gene_id_type = "SYMBOL")
p = plot_markers(marker_check)
p
ggsave(filename = paste(output_dir, "3CellMarker/marker_gene_check.database.png", sep = "/"), plot = p, width=1344,height=960,units="px", path = "./")
# 评估结果会以红色字体提示哪些marker基因在数据库中找不到对应的Ensembl名称，
# 以及哪些基因的特异性不高（标注“High overlap with XX cells”）。
# 我们可以根据评估结果优化marker基因，或者添加其他信息来辅助区分细胞类型。

# 2.4 训练分类器
# 使用marker file和cds对象训练分类器 # 这一步比较慢
sc_seurat_obj_classifier <- train_cell_classifier(cds = cds,
                                                  marker_file = marker_file,
                                                  db = org.Hs.eg.db,
                                                  cds_gene_id_type = "SYMBOL",
                                                  num_unknown = 10,
                                                  marker_file_gene_id_type = "SYMBOL",
                                                  min_observations = 50,
                                                  # cores = 8,  # windows
                                                  cores = 64, # linux
                                                  )

# 查看分类器最后选择的根节点基因，注意markerfile的基因都会在其中
feature_genes_root <- get_feature_genes(sc_seurat_obj_classifier, node="root", db= org.Hs.eg.db)
# head(feature_genes_root)

# 查看分类器中分支节点的分类基因
feature_genes_branch <- get_feature_genes(sc_seurat_obj_classifier, node = "root", db= org.Hs.eg.db, convert_ids = TRUE)
# head(feature_genes_branch)

# 3. 使用训练好的分类器预测自己的数据
pData(cds)$garnett_cluster <- pData(cds)$seurat_clusters
# 使用前面训练的pbmc_classifier来对自己的数据进行细胞分型
cds <- classify_cells(cds, sc_seurat_obj_classifier, db=org.Hs.eg.db, cluster_extend = TRUE, cds_gene_id_type = "SYMBOL")
## 将结果返回给seurat对象# 提取分类结果
cds.meta <- subset(pData(cds), select = c("cell_type","cluster_ext_type")) %>% as.data.frame()
# cds.meta <- subset(pData(cds), select = c("cell_type", "cluster_ext_type")) %>% as.data.frame()
sc_seurat_obj <- AddMetaData(sc_seurat_obj, metadata = cds.meta)

# 查看结果
p <- DimPlot(sc_seurat_obj, group.by = "cluster_ext_type", label = T, label.size = 3, repel = TRUE) + ggtitle("Classified by Garnett")
p
ggsave(filename = paste(output_dir, "4Anntotion/Garnett/Garnett.png", sep = "/"), plot = p, width=1344,height=960,units="px", path = "./")


# this is for customized.
# sc_seurat_obj@meta.data$singleR.target = sc_seurat_obj_cellType[match(sc_seurat_obj_clusters, sc_seurat_obj_cellType$ClusterID),'mona']
# sc_seurat_obj$singleR.target[which(sc_seurat_obj$clonotype_id == "clonotype42" )] <- "clon42is215"
# sc_seurat_obj$singleR.target[which(sc_seurat_obj$clonotype_id == "clonotype41" )] <- "clon41is1078"
# p5_t = DimPlot(sc_seurat_obj, reduction = "umap", label = T, group.by = 'singleR.target',)+ ggtitle("MonacoImmuneData")
# ggsave(filename=paste(output_dir,"P5_t.cluster.annotation.with.targ.png",sep="/"),plot = p5_t,width=9,height=7, path = ".")

#######################################
#######################################
write.table(data.frame("barcode"=rownames(sc_seurat_obj@meta.data),sc_seurat_obj@meta.data),
            file=paste(output_dir, 'Total.meta.data.mass.csv', sep="/"),
            sep=",", row.names=FALSE)

#########################################
############### for  CITEseq
Key(sc_seurat_obj[["ADT"]])
# p6_1 <- FeaturePlot(sc_seurat_obj, "adt_CITE-CD4", cols = c("lightgrey", "darkgreen")) + ggtitle("CD4 antibody")
# p6_2 <- FeaturePlot(sc_seurat_obj, "adt_CITE-CD8", cols = c("lightgrey", "darkgreen")) + ggtitle("CD8 antibody")
# p6_3 <- FeaturePlot(sc_seurat_obj, "adt_CITE-IgG1", cols = c("lightgrey", "darkgreen")) + ggtitle("IgG1 antibody")
# p6_4 <- FeaturePlot(sc_seurat_obj, "adt_CITE-CD3", cols = c("lightgrey", "darkgreen")) + ggtitle("CD3 antibody")
# p6_5 <- FeaturePlot(sc_seurat_obj, "CD4") + ggtitle("CD4 RNA")
# p6_6 <- FeaturePlot(sc_seurat_obj, "CD8A") + ggtitle("CD8A RNA")
# p6_7 <- FeaturePlot(sc_seurat_obj, "CD8B") + ggtitle("CD8B RNA")
# p6_8 <- FeaturePlot(sc_seurat_obj, "CD3E") + ggtitle("CD3E RNA")
library(cowplot)
# p6 = cowplot::plot_grid(p6_1,p6_2,p6_3,p6_4,p6_5,p6_6,p6_7,p6_8,ncol=4)
# p6
ggsave(filename = paste(output_dir, "P6.FeaturePlot_ADT_VS_RNA_CD4_VS_CD8.png",sep="/"), plot = p6, path = "./",width=1344,height=960,units="px",)
# question: how to make scale same.
# test:
CITE_seq_genes = c("adt_CITE-CD4","adt_CITE-CD8","adt_CITE-IgG1","adt_CITE-CD3")
RNA_seq_genes = c( "CD4","CD8A","CD8B","CD3E")
p6_1 = FeaturePlot(sc_seurat_obj, features = CITE_seq_genes, keep.scale = "all", ncol = 1)
p6_2 = FeaturePlot(sc_seurat_obj, features = RNA_seq_genes, keep.scale = "all", ncol = 1)
p6 = cowplot::plot_grid(p6_1,p6_2)
p6
ggsave(filename = paste(output_dir, "P6.FeaturePlot_ADT_VS_RNA.all.png",sep="/"), plot = p6, path = "./",width=1344,height=960,units="px",)

p6_3 = FeaturePlot(sc_seurat_obj, features = CITE_seq_genes, keep.scale = "feature", ncol = 1)
p6_4 = FeaturePlot(sc_seurat_obj, features = RNA_seq_genes, keep.scale = "feature", ncol = 1)
p6_5 = cowplot::plot_grid(p6_3, p6_4)
p6_5
ggsave(filename = paste(output_dir, "P6.FeaturePlot_ADT_VS_RNA.feature.png",sep="/"), plot = p6, path = "./",width=1344,height=960,units="px",)

# p7_1 = VlnPlot(sc_seurat_obj, "adt_CITE-CD4")
# p7_2 = VlnPlot(sc_seurat_obj, "adt_CITE-CD8")
# p7 = cowplot::plot_grid(p7_1,p7_2,ncol = 1)
# p7
p7_1 = VlnPlot(sc_seurat_obj, features = CITE_seq_genes, ncol = 1)
p7_2 = VlnPlot(sc_seurat_obj, features = RNA_seq_genes, ncol = 1)
p7 = cowplot::plot_grid(p7_1, p7_2)
p7
ggsave(filename = paste(output_dir, "P7_VlnPlot_adt.png",sep="/"), plot = p7, path = "./",width=1344,height=960,units="px",)

# we can also identify alternative protein and RNA markers for this cluster through differential
# expression
# adt_markers <- FindMarkers(sc_seurat_obj, ident.1 = 4, assay = "ADT")
# rna_markers <- FindMarkers(sc_seurat_obj, ident.1 = 1, assay = "RNA")
# adt_markers
# rna_markers
# Why do this?
# FeatureScatter(sc_seurat_obj, feature1 = "adt_CD4-TotalSeqC", feature2 = "adt_CD8-TotalSeqC")
# FeatureScatter(sc_seurat_obj, feature1 = "adt_CD4-TotalSeqC", feature2 = "rna_CD4")
# FeatureScatter(sc_seurat_obj, feature1 = "adt_CD4-TotalSeqC", feature2 = "adt_CD8-TotalSeqC")

# Ridge plots are another useful visualization.
p8_1 = RidgePlot(sc_seurat_obj, features = CITE_seq_genes, ncol = 1)
p8_2 = RidgePlot(sc_seurat_obj, features = RNA_seq_genes, ncol = 1)
p8 = plot_grid(p8_1,p8_2)
p8
ggsave(filename = paste(output_dir, "P8_1_RidgePlot_adt_VS_RNA.png",sep="/"), plot = p8, path = "./",width=1344,height=960,units="px",scaling=0.5)

# adt.markers <- FindAllMarkers(sc_seurat_obj, assay = "ADT", only.pos = TRUE)
# mean?
adt.markers
#save.image(file = paste(output_dir,"sc_seurat_obj.RData",sep="/"), version = NULL, ascii = FALSE, safe = TRUE)

# source : https://satijalab.org/seurat/archive/v3.1/multimodal_vignette.html
# source : https://satijalab.org/seurat/articles/multimodal_vignette.html


