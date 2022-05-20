setwd("/cygene2/work/P0000-Blackbird/2103-HC002/HC002004/G471_E1E2E3_scRNAseq_integration_analysis_dsy")
library(Seurat)
library(sctransform)
library(tidyverse, quietly=T)
args = commandArgs(T)

## setting for linux
# scRNAseq_path = args[1]   # scRNAseq data must contain ADT information.
# tcr_path = args[2]
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

RNA_data_E1 = data1$`Gene Expression` %>% as.matrix() %>% t() %>%  as.data.frame() %>%   rownames_to_column("barcode") %>% 
   select(c("barcode","CD4","CD8A","CD8B","CD3E","CD3D","IL7R","NKG7"))
RNA_data_E2 = data2$`Gene Expression` %>% as.matrix() %>% t() %>%  as.data.frame() %>%   rownames_to_column("barcode") %>% 
  select(c("barcode","CD4","CD8A","CD8B","CD3E","CD3D","IL7R","NKG7"))
RNA_data_E3 = data3$`Gene Expression` %>% as.matrix() %>% t() %>%  as.data.frame() %>%   rownames_to_column("barcode") %>% 
  select(c("barcode","CD4","CD8A","CD8B","CD3E","CD3D","IL7R","NKG7"))

RNA_ADT_merge_E1 = merge(RNA_data_E1,ADT_data_E1)
RNA_ADT_merge_E2 = merge(RNA_data_E2,ADT_data_E2)
RNA_ADT_merge_E3 = merge(RNA_data_E3,ADT_data_E3)

write.table(RNA_ADT_merge_E1, file=paste(output_dir,'HC002004_G471E1.RNA.ADT.csv',sep="/"), sep = ",",row.names = FALSE)
write.table(RNA_ADT_merge_E2, file=paste(output_dir,'HC002004_G471E2.RNA.ADT.csv',sep="/"), sep = ",",row.names = FALSE)
write.table(RNA_ADT_merge_E3, file=paste(output_dir,'HC002004_G471E3.RNA.ADT.csv',sep="/"), sep = ",",row.names = FALSE)

sc1_obj <- CreateSeuratObject(data1$`Gene Expression`, project = "E1")
sc2_obj <- CreateSeuratObject(data2$`Gene Expression`, project = "E2") 
sc3_obj <- CreateSeuratObject(data3$`Gene Expression`, project = "E3") 
#####################################
#############################
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
  tcr_data <- merge(tcr_data, clono[, c("clonotype_id", "cdr3s_aa","frequency","proportion")])
  # Reorder so barcodes are first column and set them as rownames.
  tcr_data <- tcr_data %>% tibble::column_to_rownames("barcode")
  
  # Add to the Seurat object's metadata.
  clono_seurat <- Seurat::AddMetaData(object=seurat_obj, metadata=tcr_data)
  return(clono_seurat)
}
sc1_obj = add_clonotype(TCR_path1, sc1_obj)
sc2_obj = add_clonotype(TCR_path2, sc2_obj)
sc3_obj = add_clonotype(TCR_path3, sc3_obj)

# 数据整合：

########################################
#########################################
#################  integration
################# 
sc.list <- list("E1" = sc1_obj, "E2" = sc2_obj,"E3" = sc3_obj)
sc.list <- lapply(X=sc.list, FUN = function(x){
  x <- PercentageFeatureSet(x, pattern = "^MT", col.name = "percent.mt")
  x <- SCTransform(x, method = "glmGamPoi", vars.to.regress = "percent.mt", verbose = FALSE)
})
features <- SelectIntegrationFeatures(object.list = sc.list, nfeatures=3000)
sc.list <- PrepSCTIntegration(object.list = sc.list, anchor.features = features)

immune.anchors <- FindIntegrationAnchors(object.list = sc.list, normalization.method = "SCT", anchor.features = features) 
immune.combined.sct <- IntegrateData(anchorset = immune.anchors, normalization.method = "SCT") 
immune.combined.sct <- RunPCA(immune.combined.sct) %>% RunUMAP(dims = 1:30) %>% 
  FindNeighbors(dims=1:30) %>% FindClusters(resolution=0.5)

############################################ UMAP 
DimPlot(immune.combined.sct, label=T) 
p = FeaturePlot(immune.combined.sct, features = c("CD4","CD8A","FOXP3","CCR7","IL26","EOMES","KLRD1","IFNG","PDCD1")) +
  theme(
    panel.background = element_rect(fill = 'transparent'), # bg of the panel
    plot.background = element_rect(fill = 'transparent', color = NA), # bg of the plot
    panel.grid.major = element_blank(), # get rid of major grid
    panel.grid.minor = element_blank(), # get rid of minor grid
    legend.background = element_rect(fill = 'transparent'), # get rid of legend bg
    legend.box.background = element_rect(fill = 'transparent'), # get rid of legend panel bg
    legend.key = element_rect(fill = "transparent", colour = NA), # get rid of key legend fill, and of the surrounding
    axis.line = element_line(colour = "black") # adding a black line for x and y axis
  )
ggsave(filename = "myplot.png", p, device = "png", width =1200, height = 900,units = "px",dpi = 72, bg="transparent")
ggsave(filename = "myplot.png", p, device = "png", width =120, height = 90,  units = "cm",dpi = 72, bg="transparent")
#########################################




# merged_obj = merge(sc1_obj, y=c(sc2_obj, sc3_obj),add.cell.ids = c("E1","E2","E3"), project = "G471")
# merged_obj
# merged_obj = PercentageFeatureSet(merged_obj, pattern = "^MT-", col.name = "percent.mt") %>% 
  # ScaleData() %>% FindVariableFeatures() %>% RunPCA() %>% RunUMAP(dims=1:30) %>% 
  # FindNeighbors(dims = 1:30) %>% FindClusters(resolution = 0.5)
# before_integration_umap_sample = DimPlot(merged_obj, group.by = "orig.ident")

# doubletsfinder
#########################################
# find doublets
## install.packages('remotes')
## remotes::install_github('chris-mcginnis-ucsf/DoubletFinder')
library(DoubletFinder)
sc_seurat_obj = immune.combined.sct
# doublets_detecte <- function(sc_seurat_obj){
# Do not apply DoubletFinder to aggregated scRNA-seq data representing multiple distinct samples
## -------------- pK Identification (no ground-truth) --------------------------------------
sweep.res.list <- DoubletFinder::paramSweep_v3(sc_seurat_obj, PCs = 1:30, sct = TRUE, num.cores = 64)
sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)
bcmvn <- find.pK(sweep.stats)
bcmvn
mpK <- as.numeric(as.vector(bcmvn$pK[which.max(bcmvn$BCmetric)]))
# mpK ?
## ----------------- DoubletFinder:Homotypic Doublet Proportion Estimate ------------------------
annotations <- sc_seurat_obj@meta.data$seurat_clusters
homotypic.prop <- modelHomotypic(annotations)
nExp_poi <- round(0.075*nrow(sc_seurat_obj@meta.data)) ## Assuming 7.5% doublet formation rate
# DoubletRate = 0.039  # 5000 Cells correspond to doublets rate yes 3.9%
# nExp_poi <- round(DoubletRate*ncol(sc_seurat_obj))
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

## ------------DoubletFinder:Run DoubletFinder with varying classification stringencies ---------
sc_seurat_obj <- doubletFinder_v3(sc_seurat_obj, PCs = 1:30, pN = 0.25, pK = mpK,
                                  nExp = nExp_poi, reuse.pANN = FALSE, sct = TRUE)
sc_seurat_obj <- doubletFinder_v3(sc_seurat_obj, PCs = 1:30, pN = 0.25, pK = mpK,
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
p4

########################## finish Doubletsfinder.
##########################################################
#########################################################
# QC
VlnPlot(immune.combined.sct,features = c("nFeature_RNA", "nCount_RNA","percent.mt"), group.by = "orig.ident")
# immune.combined.sct = subset(immune.combined.sct)
# Visualization
p1 <- DimPlot(immune.combined.sct, reduction = "umap",group.by = "orig.ident") +ggtitle("")
p2 <- DimPlot(immune.combined.sct, reduction = "umap", label = TRUE, repel = TRUE)

p1# 图没保存
p2
T_cell_marker = c("CD3E","CD4","IL7R","CD8A","CD8B","NKG7")
FeaturePlot(immune.combined.sct, features = T_cell_marker, label = T)
DotPlot(immune.combined.sct, features = T_cell_marker,)
RidgePlot(immune.combined.sct, features = T_cell_marker)
GetAssay(immune.combined.sct)
B_cell_marker = c("CD19","CD79A")

RidgePlot(immune.combined.sct, features = B_cell_marker)
# DefaultAssay(immune.combined.sct) <- "RNA"
# DefaultAssay(immune.combined.sct) <- "integrated"
DefaultAssay(immune.combined.sct) <- "SCT"
GetAssay(immune.combined.sct)
CD4_barcode = WhichCells(immune.combined.sct, expression = CD4>0, slot = "data",)
length(CD4_barcode)

CD8_barcode = WhichCells(immune.combined.sct, expression = CD8A>0 | CD8B>0)
length(CD8_barcode)
immune.combined.sct$manual_label_cells <- "Others"
immune.combined.sct$manual_label_cells[rownames(immune.combined.sct@meta.data) %in% CD4_barcode] <- "CD4"
immune.combined.sct$manual_label_cells[rownames(immune.combined.sct@meta.data) %in% CD8_barcode] <- "CD8"

DimPlot(immune.combined.sct, reduction = "umap", group.by = "manual_label_cells", cols = c("red","blue","grey80"), pt.size = 0.3) + ggtitle("")


(immune.combined.sct)

# prop.table(table(Idents(merged_obj)))
# write.table(table(Idents(merged_obj), merged_obj$orig.ident), file = "mytest.csv", sep="\t", row.names = T)
# write.table(prop.table(table(Idents(merged_obj))),file = "mytest.prop.csv", sep = "\t",row.names = T)


# ################# after QC
after_obj <- subset(immune.combined.sct, subset = nFeature_RNA > 100 & percent.mt < 10)
after_obj
CD4_barcode = WhichCells(after_obj, expression = CD4>0, slot = "data",)
length(CD4_barcode)
CD8_barcode = WhichCells(after_obj, expression = CD8A>0 | CD8B>0)
length(CD8_barcode)

RidgePlot(immune.combined.sct, features = T_cell_marker, log = FALSE)
RidgePlot(immune.combined.sct, features = T_cell_marker, log = TRUE)
bcell_marker = c("CD19","CD79A")
RidgePlot(immune.combined.sct, features = bcell_marker, log=T)
table(Idents(immune.combined.sct))
table(Idents(immune.combined.sct), immune.combined.sct$orig.ident)

ggsave(path = output_dir, filename = "integration.UMAP.png",plot = p3, width = 11, height = 7,bg = "transparent",device = "png")
# To visualize the two conditions side-by-side, we can use the split.by argument to show each condition colored by cluster.
# DimPlot(immune.combined.sct, reduction = "umap", split.by = "orig.ident")

T_cells_clusters = c(1,4,6,8,12,14)

T_cells_barcodes = WhichCells(immune.combined.sct, idents = T_cells_clusters)
head(T_cells_barcodes)


immune.combined.sct$after_cell_type <- "Non-T Cell"
immune.combined.sct$after_cell_type[rownames(immune.combined.sct@meta.data) %in% T_cells_barcodes] <- "T-cell"
head(immune.combined.sct@meta.data)
write.table(table(Idents(immune.combined.sct), immune.combined.sct$after_cell_type), file = "Tcell_number_of_cell_type_after_qc.xls", sep = "\t", row.names = T)

DimPlot(immune.combined.sct, group.by = "after_cell_type")+ ggtitle("cell type (after QC)")
table(Idents(immune.combined.sct), immune.combined.sct$after_cell_type)

####### 给cluster重命名：
# 现在已知T 细胞的cluster为：（1,4,6,8,12,14）
# B 细胞的cluster为15.
immune.combined.sct <- RenameIdents(immune.combined.sct, 
                                    "1" = "T cells","4"="T cells",
                                    "6"="T cells","8"="T cells",
                                    "12"="T cells","14"="T cells",
                                    "15"="B cells")
DimPlot(immune.combined.sct,label = T)






######
FeaturePlot(immune.combined.sct, features = c("CD4"))
  theme(
    axis.text  = element_blank(), # keduwenzi
    axis.ticks = element_blank(), # kedu
    axis.line  = element_blank(),
    panel.border = element_rect(fill=NA, color = "black", size=1, linetype = "solid")
  ) +xlab("")+ylab("")+
  NoLegend()
#  
other_celltype = c("PTPRC","CD14","CD68","CD163","ITGAX","ITGAM","CD33",
                   "CD19","CD79A","NCAM1","FCGR3A", "CD3E","CD4","IL7R",
                   "CD8A","NKG7","CD4","IL2RA","FOXP3","PECAM1","CD34")
tumor_marker = c("AFP","GPC3")
LCSC_marker = c("ALDH1A1","EPCAM","KRT19","ANPEP","CD24","CD44","CD47","THY1","PROM1")

umap_data = Embeddings(immune.combined.sct, reduction = "umap")
# 
# DefaultAssay(immune.combined.sct) <- ""
FeaturePlot(immune.combined.sct, features = other_celltype)
FeaturePlot(immune.combined.sct,features = tumor_marker)
FeaturePlot(immune.combined.sct, features = LCSC_marker)
DefaultAssay(immune.combined.sct) <- "RNA"
RidgePlot(immune.combined.sct,features=tumor_marker,)

library(SingleR)
library(celldex)
DefaultAssay(immune.combined.sct) <- "SCT"
sc_seurat_obj_SingleR = GetAssayData(immune.combined.sct, slot="data")
sc_seurat_obj_clusters = immune.combined.sct@meta.data$seurat_clusters

hpca = celldex::HumanPrimaryCellAtlasData()
blueprint_encode = celldex::BlueprintEncodeData()
dice = celldex::DatabaseImmuneCellExpressionData()
mona = celldex::MonacoImmuneData()
novershter = celldex::NovershternHematopoieticData()

sc_seurat_obj = immune.combined.sct
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
p5_1
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
as.character(immune.combined.sct@meta.data$seurat_clusters)
scCATCH_obj <- scCATCH::createscCATCH(immune.combined.sct[["RNA"]]@data, cluster = as.character(immune.combined.sct@meta.data$seurat_clusters))
scCATCH_obj <- scCATCH::createscCATCH(immune.combined.sct[["SCT"]]@data, cluster = as.character(immune.combined.sct@meta.data$seurat_clusters))
tumor_cluster = c(0,2,3,5,7,9,10,11,13)
clu_markers <- scCATCH::findmarkergene(scCATCH_obj, species = "Human", 
                                        cluster = tumor_cluster, cancer = "Hepatocellular Cancer", 
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
FindMarkers(immune.combined.sct, ident.1 = "1",only.pos = T, min.pct = 0.25, logfc.threshold = 0.25)
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
expression_matrix = GetAssayData(immune.combined.sct, assay="RNA", slot = 'counts')
cell_metadata <- immune.combined.sct@meta.data
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
marker_file = "/cygene2/work/P0000-Blackbird/2103-HC002/HC002004/G471_E1E2E3_scRNAseq_integration_analysis_dsy/HC_markers.txt"
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
immune.combined.sct <- AddMetaData(immune.combined.sct, metadata = cds.meta)

# 查看结果
p <- DimPlot(immune.combined.sct, group.by = "cluster_ext_type", label = T, label.size = 3, repel = TRUE) + ggtitle("Classified by Garnett")
p
ggsave(filename = paste(output_dir, "4Annotation/Garnett/Garnett.png", sep = "/"), plot = p, width=1344,height=960,units="px", path = "./")


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
#
#########################################

#########################################
############### for  CITEseq
#################################################
sc1_obj[['ADT']] <- CreateAssayObject(counts = data1$`Antibody Capture`)
sc2_obj[['ADT']] <- CreateAssayObject(counts = data2$`Antibody Capture`)
sc3_obj[['ADT']] <- CreateAssayObject(counts = data3$`Antibody Capture`)

rownames(sc1_obj[["ADT"]])
## [1] "CD8-CITE" "CD4-CITE" "CD3-CITE" "IgG1-CITE"

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
#####################################################################################################
# troubleshooting
sc1_obj <- PercentageFeatureSet(sc1_obj, pattern = "^MT-", col.name = "percent.mt") %>%
  SCTransform(method = "glmGamPoi", vars.to.regress = "percent.mt", verbose = FALSE) %>%
  RunPCA(npcs=30) %>%
  FindNeighbors(dims = 1:10) %>%
  RunUMAP(dims = 1:10) %>%
  FindClusters(resolution=0.5)

DimPlot(sc1_obj,label = T)
DefaultAssay(sc1_obj) # is SCT
sc1_obj

FeaturePlot(sc1_obj, features = c("CD4","CD8A","CD8B"))

TCRseq_barcode = rownames(sc1_obj@meta.data[which(sc1_obj@meta.data$clonotype_id != "NA"),])
length(TCRseq_barcode)

# DimPlot(sc1_obj, reduction = "umap", cells.highlight = TCRseq_barcode, pt.size = 0.3)+
#   scale_color_manual(labels = c("Others","TCR cells"),values = c("grey50","blue")) #  + labs(color = "legend title")

sc1_obj$TCR_cells <- "Others"
sc1_obj$TCR_cells[rownames(sc1_obj@meta.data) %in% TCRseq_barcode] <- "TCR Cells"
DimPlot(sc1_obj, reduction = "umap", group.by = "TCR_cells", cols = c("grey50","blue"),pt.size = 0.2) + ggtitle("")

#########################################
############### for  CITEseq
#################################################
sc1_obj[['ADT']] <- CreateAssayObject(counts = data1$`Antibody Capture`)

rownames(sc1_obj[["ADT"]])
## [1] "CD8-CITE" "CD4-CITE" "CD3-CITE" "IgG1-CITE"

DefaultAssay(sc1_obj) <- "ADT"

# CITE_seq_genes = c("CD8-CITE", "CD4-CITE", "IgG1-CITE")
CITE_seq_genes = rownames(sc1_obj[["ADT"]])
CITE_seq_genes
FeaturePlot(sc1_obj, features = CITE_seq_genes, reduction = "umap", keep.scale = "all")
FeaturePlot(sc1_obj, features = CITE_seq_genes, reduction = "umap", keep.scale = "feature")

