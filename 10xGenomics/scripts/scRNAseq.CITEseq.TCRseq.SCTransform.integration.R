# install.packages("Seurat")
# install.packages("tidyr")
# install.packages("tidyverse")
library(Seurat)
library(sctransform)
library(tidyverse, quietly=T)
args = commandArgs(T)


## setting for linux
# scRNAseq_path = args[1]   # scRNAseq data must contain ADT information.
# tcr_path = args[2]


# create a output directory in current path.
dir_list <- c("1QC",
              "2Cluster",
              "3CellMarker",
              "4Annotation",
              "4Annotation/SingleR",
              "4Annotation/Garnett",
              "4Annotation/scCATCH",
              "5TrajectoryAnalysis",
              "6RNAvelocity",
              "7CellPhone",
              "8scCNV")
for (each in dir_list){
  if(!dir.exists(each)){
    dir.create(each)
  }
}


output_dir = "."
# read scRNAseq folder:
scRNAseq_path1 = "/cygene2/work/P0000-Blackbird/2103-HC002/HC002004/HC24_G471E1L2_scRNAseq_G471E1L3_CITEseq/outs"
scRNAseq_path2 = "/cygene2/work/P0000-Blackbird/2103-HC002/HC002004/HC24_G471E2L2_scRNAseq_G471E2L3_CITEseq/outs"
scRNAseq_path3 = "/cygene2/work/P0000-Blackbird/2103-HC002/HC002004/HC24_G471E3L2_scRNAseq_G471E3L3_CITEseq/outs"
data1 = Seurat::Read10X(data.dir = paste(scRNAseq_path1, '/filtered_feature_bc_matrix', sep = "/"))
data2 = Seurat::Read10X(data.dir = paste(scRNAseq_path2, '/filtered_feature_bc_matrix', sep = "/"))
data3 = Seurat::Read10X(data.dir = paste(scRNAseq_path3, '/filtered_feature_bc_matrix', sep = "/"))

# 10X data contains more than one type and is being returned as a list containing matrices of each type.

ADT_data_E1 = data1$`Antibody Capture` %>% as.matrix() %>% t() %>% as.data.frame() %>% rownames_to_column("barcode")
ADT_data_E2 = data2$`Antibody Capture` %>% as.matrix() %>% t() %>% as.data.frame() %>% rownames_to_column("barcode")
ADT_data_E3 = data3$`Antibody Capture` %>% as.matrix() %>% t() %>% as.data.frame() %>% rownames_to_column("barcode")
write.table(ADT_data_E1, file=paste(output_dir,'HC24_G471E1L3_ADT.t.csv',sep="/"), sep = ",",row.names = FALSE)
write.table(ADT_data_E2, file=paste(output_dir,'HC24_G471E2L3_ADT.t.csv',sep="/"), sep = ",",row.names = FALSE)
write.table(ADT_data_E3, file=paste(output_dir,'HC24_G471E3L3_ADT.t.csv',sep="/"), sep = ",",row.names = FALSE)

RNA_data_E1 = data1$`Gene Expression` %>% as.matrix() %>% t() %>%  as.data.frame() %>%   rownames_to_column("barcode") %>% 
   select(c("barcode","CD4","CD8A","CD8B","CD3E","CD3D"))
RNA_data_E2 = data2$`Gene Expression` %>% as.matrix() %>% t() %>%  as.data.frame() %>%   rownames_to_column("barcode") %>% 
  select(c("barcode","CD4","CD8A","CD8B","CD3E","CD3D"))
RNA_data_E3 = data3$`Gene Expression` %>% as.matrix() %>% t() %>%  as.data.frame() %>%   rownames_to_column("barcode") %>% 
  select(c("barcode","CD4","CD8A","CD8B","CD3E","CD3D"))
write.table(RNA_data_E1, file=paste(output_dir,'HC24_G471E1L2.RNA.t.csv',sep="/"), sep = ",", row.names = FALSE)
write.table(RNA_data_E2, file=paste(output_dir,'HC24_G471E2L2.RNA.t.csv',sep="/"), sep = ",", row.names = FALSE)
write.table(RNA_data_E3, file=paste(output_dir,'HC24_G471E3L2.RNA.t.csv',sep="/"), sep = ",", row.names = FALSE)

sc1_obj <- CreateSeuratObject(data1$`Gene Expression`, project = "E1")
sc2_obj <- CreateSeuratObject(data2$`Gene Expression`, project = "E2") 
sc3_obj <- CreateSeuratObject(data3$`Gene Expression`, project = "E3") 

sc1_obj[['ADT']] <- CreateAssayObject(counts = data1$`Antibody Capture`)
sc2_obj[['ADT']] <- CreateAssayObject(counts = data2$`Antibody Capture`)
sc3_obj[['ADT']] <- CreateAssayObject(counts = data3$`Antibody Capture`)

rownames(sc1_obj[["ADT"]])
## [1] "CD8-CITE" "CD4-CITE" "CD3-CITE" "IgG1-CITE"

################
sc1_obj = PercentageFeatureSet(sc1_obj, pattern = "^MT-", col.name = "percent.mt")
sc2_obj = PercentageFeatureSet(sc2_obj, pattern = "^MT-", col.name = "percent.mt")
sc3_obj = PercentageFeatureSet(sc3_obj, pattern = "^MT-", col.name = "percent.mt")

sc1_obj = SCTransform(sc1_obj,method = "glmGamPoi", vars.to.regress = "percent.mt",verbose=FALSE)
sc2_obj = SCTransform(sc2_obj,method = "glmGamPoi", vars.to.regress = "percent.mt",verbose=FALSE)
sc3_obj = SCTransform(sc3_obj,method = "glmGamPoi", vars.to.regress = "percent.mt",verbose=FALSE)


#########################################
#########################################
# for TCR seq
TCR_path1 = "/cygene2/work/P0000-Blackbird/2103-HC002/HC002004/HC24_G471E1L1_TCRseq_IMGT"
TCR_path2 = "/cygene2/work/P0000-Blackbird/2103-HC002/HC002004/HC24_G471E2L1_TCRseq_IMGT"
TCR_path3 = "/cygene2/work/P0000-Blackbird/2103-HC002/HC002004/HC24_G471E3L1_TCRseq_IMGT"

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
sc1_obj = add_clonotype(TCR_path1, sc1_obj)
sc2_obj = add_clonotype(TCR_path2, sc2_obj)
sc3_obj = add_clonotype(TCR_path3, sc3_obj)


#########################################
#########################################


##### integration
sc.list <- list('E1'=sc1_obj, "E2" = sc2_obj, "E3" = sc3_obj)
features <- SelectIntegrationFeatures(object.list = sc.list, nfeatures=3000)
sc.list <- PrepSCTIntegration(object.list = sc.list, anchor.features = features)

immune.anchors <- FindIntegrationAnchors(object.list = sc.list, normalization.method = "SCT", anchor.features = features) # In CheckDuplicateCellNames(object.list = object.list) :
# Some cell names are duplicated across objects provided. Renaming to enforce unique cell names.
immune.combined.sct <- IntegrateData(anchorset = immune.anchors, normalization.method = "SCT") 
immune.combined.sct <- RunPCA(immune.combined.sct) %>% RunUMAP(dims = 1:30) %>% 
  FindNeighbors(dims=1:30) %>% FindClusters(resolution=0.5)

# Visualization
p1 <- DimPlot(immune.combined.sct, reduction = "umap",group.by = "orig.ident")
p2 <- DimPlot(immune.combined.sct, reduction = "umap", label = TRUE, repel = TRUE)
p1 + p2
ggsave()
# To visualize the two conditions side-by-side, we can use the split.by argument to show each condition colored by cluster.
DimPlot(immune.combined.sct, reduction = "umap", split.by = "orig.ident")
ggsave()

barcode = rownames(immune.combined.sct@meta.data[which(immune.combined.sct@meta.data$clonotype_id != "NA"),])
DimPlot(immune.combined.sct, reduction = "umap", cells.highlight = barcode,split.by = "orig.ident")
ggsave()
# change unselected, and group_1
######




library(SingleR)
library(celldex)
DefaultAssay(immune.combined.sct) <- "RNA"
sc_seurat_obj_SingleR = GetAssayData(immune.combined.sct, slot="data")
sc_seurat_obj_clusters = immune.combined.sct@meta.data$seurat_clusters

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

sc_seurat_obj_cellType = data.frame(ClusterID=levels(sc_seurat_obj_clusters),
                                    hpca = pred.hpca$labels,
                                    blue  = pred.blue$labels,
                                    Dice = pred.dice$labels,
                                    mona = pred.mona$labels,
                                    novershter = pred.novershter$labels
)
sc_seurat_obj_cellType
write.table(sc_seurat_obj_cellType, file = paste(output_dir,"P5.SingleR_cell_type.csv",sep="/"),sep = ",", row.names = FALSE, quote=FALSE)
immune.combined.sct@meta.data$singleR.hpca = sc_seurat_obj_cellType[match(sc_seurat_obj_clusters, sc_seurat_obj_cellType$ClusterID),'hpca']
p5_0 = DimPlot(immune.combined.sct, reduction = "umap", label = T, group.by = 'singleR.hpca',)+ggtitle("HumanPrimaryCellAtlasData")
p5_0
immune.combined.sct@meta.data$singleR.blue = sc_seurat_obj_cellType[match(sc_seurat_obj_clusters, sc_seurat_obj_cellType$ClusterID),'blue']
p5_1 = DimPlot(immune.combined.sct, reduction = "umap", label = T, group.by = 'singleR.blue',)+ggtitle("BlueprintEncodeData")
immune.combined.sct@meta.data$singleR.dice = sc_seurat_obj_cellType[match(sc_seurat_obj_clusters, sc_seurat_obj_cellType$ClusterID),'Dice']
p5_2 = DimPlot(immune.combined.sct, reduction = "umap", label = T, group.by = 'singleR.dice',)+ggtitle("DatabaseImmuneCellExpressionData")
immune.combined.sct@meta.data$singleR.mona = sc_seurat_obj_cellType[match(sc_seurat_obj_clusters, sc_seurat_obj_cellType$ClusterID),'mona']
p5_3 = DimPlot(immune.combined.sct, reduction = "umap", label = T, group.by = 'singleR.mona') + ggtitle("MonacoImmuneData")
immune.combined.sct@meta.data$singleR.nove = sc_seurat_obj_cellType[match(sc_seurat_obj_clusters, sc_seurat_obj_cellType$ClusterID),'novershter']
p5_4 = DimPlot(immune.combined.sct, reduction = "umap", label = T, group.by = 'singleR.nove',)+ ggtitle("NovershternHematopoieticData")
ggsave(filename=paste(output_dir,"4Annotation/SingleR/P5_0.cluster.annotation.with.hpca.png",sep="/"),width=1344,height=960,units="px",plot = p5_0, path = "./")
ggsave(filename=paste(output_dir,"4Annotation/SingleR/P5_1.cluster.annotation.with.blue.png",sep="/"),width=1344,height=960,units="px",plot = p5_1, path = "./")
ggsave(filename=paste(output_dir,"4Annotation/SingleR/P5_2.cluster.annotation.with.dice.png",sep="/"),width=1344,height=960,units="px",plot = p5_2, path = "./")
ggsave(filename=paste(output_dir,"4Annotation/SingleR/P5_3.cluster.annotation.with.mona.png",sep="/"),width=1344,height=960,units="px",plot = p5_3, path = "./")
ggsave(filename=paste(output_dir,"4Annotation/SingleR/P5_4.cluster.annotation.with.nove.png",sep="/"),width=1344,height=960,units="px",plot = p5_4, path = "./")
library(cowplot)
p5 = cowplot::plot_grid(p5_0,p5_1,p5_2,p5_3,p5_4, ncol = 2)
p5
# p5 = cowplot::plot_grid(p5_0,p5_1,p5_3,p5_4, ncol = 2)
ggsave(filename=paste(output_dir,"4Annotation/SingleR/P5_cluster.annotation.with.5.databases.png",sep="/"), plot = p5, width=1344,height=960,units="px", path = "./" )

########################################
########################################
######## scCATCH
library(scCATCH)
scCATCH_obj <- scCATCH::createscCATCH(immune.combined.sct[["RNA"]]@data, cluster = as.character(immune.combined.sct@meta.data$seurat_clusters))
clu_markers <- scCATCH::findmarkergene(scCATCH_obj, species = "Human", 
                                        cluster = "All", cancer = "Hepatocellular Cancer", 
                                        marker = cellmatch,
                                        tissue = c("Liver"), 
                                        use_method = "1", cell_min_pct = 0.25,
                                        logfc = 0.25, pvalue = 0.25,
                                        verbose = TRUE)
cellmatch.table = cellmatch
write.csv(cellmatch.table, file = paste(output_dir,  "scCATCH.CellMatch.database.csv", sep = "/"),row.names = FALSE, quote = FALSE)
clu_ann <- scCATCH::findcelltype(clu_markers)

scCATCH_ann_data = merge(clu_ann@meta, clu_ann@celltype,by = "cluster") %>% tibble::column_to_rownames("cell")
immune.combined.sct = AddMetaData(immune.combined.sct, scCATCH_ann_data)

p5_6 = DimPlot(immune.combined.sct, reduction = "umap", label = T, group.by="cell_type", pt.size = 0.5) + ggtitle("scCATCH Peripheral blood")
p5_6
ggsave(filename=paste(output_dir, "4Annotation/scCATCH/P5_6.cluster.annotation.with.scCATCH.cellmatch.png", sep="/"), plot = p5_6, width=1344,height=960,units="px",path = "./" )

########################################
########################################

# find markers for every cluster compared to all remaining cells, report only the positive
# ones
immune.combined.sct.markers <- FindAllMarkers(immune.combined.sct, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
top5_marker_genes = immune.combined.sct.markers %>%
  group_by(cluster) %>%
  slice_max(n = 5, order_by = avg_log2FC)
write.csv(top5_marker_genes, file = paste(output_dir, "top5_marker_genes.csv",sep = "/"), row.names = FALSE, quote = F)


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
#################################################
DefaultAssay(immune.combined.sct) <- "ADT"
# VariableFeatures(immune.combined.sct) = rownames(immune.combined.sct[["ADT"]])
# immune.combined.sct <- NormalizeData(immune.combined.sct, normalization.method = "CLR", margin=2) %>% 
#   ScaleData() %>% RunPCA() %>% RunUMAP()
# immune.combined.sct <- ScaleData(immune.combined.sct) 
# immune.combined.sct <- RunPCA(immune.combined.sct)
# 
# immune.combined.sct <- FindMultiModalNeighbors(
#   immune.combined.sct, reduction.list = list("pca","apca"),
#   dims.list = list(1:30,1:4),modality.weight.name = "RNA.weight"
# )
# 
# immune.combined.sct <- RunUMAP(immune.combined.sct,nn.name = "weighted.nn", reduction.name = "wnn.umap", reduction.key = "wnnUMAP_")
CITE_seq_genes = rownames(immune.combined.sct[["ADT"]])
CITE_seq_genes
FeaturePlot(immune.combined.sct, features = CITE_seq_genes, reduction = "umap", keep.scale = "all", ncol = 1, split.by="orig.ident")
FeaturePlot(immune.combined.sct, features = CITE_seq_genes, reduction = "umap", keep.scale = "feature", ncol = 1, split.by="orig.ident")
FeaturePlot(immune.combined.sct, features = CITE_seq_genes, reduction = "umap", keep.scale = NULL, ncol = 1, split.by="orig.ident")

VlnPlot(immune.combined.sct, assay = "ADT", features = c("CD8-CITE","IgG1-CITE","CD3-CITE","CD4-CITE"),log=T)
# ggsave(filename=paste(output_dir,"1QC/P1.VlnPlot.feature.png",sep="/"), plot=p1, width=1344,height=960, units="px", path=".")

# FeatureScatter is typically used to visualize feature-feature relationships, but can be used
# for anything calculated by the object, i.e. columns in object metadata, PC scores etc.

# sc_seurat_obj <- subset(sc_seurat_obj, subset = nFeature_RNA > 20 & nFeature_RNA < 3000 & percent.mt < 25)
# Identify conserved cell type markers
DefaultAssay(immune.combined.sct) <- "RNA"
FeaturePlot(immune.combined.sct, features = c("CD4", "CD8A", "CD8B","MS4A1"), min.cutoff = "q9", keep.scale = "feature")

VlnPlot(immune.combined.sct, assay = "RNA",group.by = "orig.ident",  
        features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), log=F, flip=FALSE)
FeatureScatter(immune.combined.sct, feature1 = "nCount_RNA", feature2 = "percent.mt",group.by = "orig.ident")
FeatureScatter(immune.combined.sct, feature1 = "nCount_RNA", feature2 = "nFeature_RNA", group.by = "orig.ident")

# ##############################
# annotation

Idents(immune.combined) <- factor(Idents(immune.combined), levels = c("HSPC", "Mono/Mk Doublets",
                                                                      "pDC", "Eryth", "Mk", "DC", "CD14 Mono", "CD16 Mono", "B Activated", "B", "CD8 T", "NK", "T activated",
                                                                      "CD4 Naive T", "CD4 Memory T"))
markers.to.plot <- c("CD3D", "CREM", "HSPH1", "SELL", "GIMAP5", "CACYBP", "GNLY", "NKG7", "CCL5",
                     "CD8A", "MS4A1", "CD79A", "MIR155HG", "NME1", "FCGR3A", "VMO1", "CCL2", "S100A9", "HLA-DQA1",
                     "GPR183", "PPBP", "GNG11", "HBA2", "HBB", "TSPAN13", "IL3RA", "IGJ", "PRSS57")
DotPlot(immune.combined, features = markers.to.plot, cols = c("blue", "red"), dot.scale = 8, split.by = "stim") +
  RotatedAxis()




p3 <- DimPlot(sc_seurat_obj, reduction = "umap")
p3_2 <- DimPlot(sc_seurat_obj, reduction = "tsne")
p3
p3_2
ggsave(filename = paste(output_dir,"/2Cluster/P3.raw.umap.png",sep="/"),plot = p3, width=1344,height=960,units="px",path = "." )
ggsave(filename = paste(output_dir,"/2Cluster/P3.raw.tsne.png",sep="/"),plot = p3_2, width=1344,height=960,units="px",path = "." )


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
ggsave(filename = paste(output_dir, "P8_1_RidgePlot_adt_VS_RNA.png",sep="/"), plot = p8, path = "./",width=1344,height=960,units="px",)

# adt.markers <- FindAllMarkers(sc_seurat_obj, assay = "ADT", only.pos = TRUE)
# mean?
adt.markers
#save.image(file = paste(output_dir,"sc_seurat_obj.RData",sep="/"), version = NULL, ascii = FALSE, safe = TRUE)

# source : https://satijalab.org/seurat/archive/v3.1/multimodal_vignette.html
# source : https://satijalab.org/seurat/articles/multimodal_vignette.html


