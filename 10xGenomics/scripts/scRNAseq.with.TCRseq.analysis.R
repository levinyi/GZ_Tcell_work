# install.packages("Seurat")
# install.packages('tidyr')
# install.packages("tidyverse")
library(Seurat)
library(sctransform)
library(tidyverse, quietly=T)
library(cowplot)
# args = commandArgs(T)
setwd("/cygene2/work/P0000-Blackbird/2103-HC002/HC002009")

## setting for linux
# scRNAseq_path = args[1]
# tcr_path = args[2]

scRNAseq_path = "/cygene2/work/P0000-Blackbird/2103-HC002/HC002009/G480E1_GEX/outs/"
tcr_path = "/cygene2/work/P0000-Blackbird/2103-HC002/HC002009/G480E1_VDJ/outs/"

# create a output directory in current path.
dir.create("scRNAseq_with_TCR_analysis")
output_dir = "scRNAseq_with_TCR_analysis"

# read scRNAseq folder:
data = Seurat::Read10X(data.dir = paste(scRNAseq_path,'filtered_feature_bc_matrix',sep = "/"))

sc_obj <- Seurat::CreateSeuratObject(data, project = "scRNAseq_object")
# sc_obj # 33538 features across 5533 samples within 1 assay 
print(sc_obj)
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
  tcr_data <- merge(tcr_data, clono[, c("clonotype_id", "cdr3s_aa","frequency","proportion")])
  # Reorder so barcodes are first column and set them as rownames.
  tcr_data <- tcr_data %>% tibble::column_to_rownames("barcode")

  # Add to the Seurat object's metadata.
  clono_seurat <- Seurat::AddMetaData(object=sc_obj, metadata=tcr_data)
  return(clono_seurat)
}
sc_obj = add_clonotype(tcr_path, sc_obj)


################
sc_obj[["percent.mt"]] <- PercentageFeatureSet(sc_obj, pattern = "^MT-")
p1 = VlnPlot(sc_obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
ggsave(filename=paste(output_dir,"P1.VlnPlot.feature.png",sep="/"), plot=p1, width=12,height=6, path=".")
# FeatureScatter is typically used to visualize feature-feature relationships, but can be used
# for anything calculated by the object, i.e. columns in object metadata, PC scores etc.

plot1 <- FeatureScatter(sc_obj, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(sc_obj, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
p2 = CombinePlots(plots = list(plot1, plot2))
ggsave(filename=paste(output_dir, "P2.combinePlots.FeatureScatter.png", sep="/"), plot=p2,width=12,height=6, path=".")

sc_obj <- subset(sc_obj, subset = nFeature_RNA > 20 & nFeature_RNA < 3000 & percent.mt < 25)
# sc_obj # 33538 features across 3549 samples within 1 assay 
# standard seurat workflow
sc_obj <- SCTransform(sc_obj, method = "glmGamPoi", vars.to.regress = "percent.mt", verbose = FALSE) %>%
  RunPCA(npcs=30) %>%
  FindNeighbors(dims = 1:10) %>%
  RunUMAP(dims = 1:10) %>%
  FindClusters(resolution=0.5) 

p3 <- DimPlot(sc_obj, reduction = "umap")
p3
ggsave(filename = paste(output_dir, "P3.raw.umap.png", sep="/"), plot = p3, width=9, height=7, path = "." )


FeaturePlot(sc_obj, features = c("CD4","CD8A","CD8B"))
###########################################################
### for TCR data
TCRseq_barcode = rownames(sc_obj@meta.data[which(sc_obj@meta.data$clonotype_id != "NA"),])
length(TCRseq_barcode)

# DimPlot(sc_obj, reduction = "umap", cells.highlight = TCRseq_barcode, pt.size = 0.3)+
#   scale_color_manual(labels = c("Others","TCR cells"),values = c("grey50","blue")) #  + labs(color = "legend title")
# or
sc_obj$TCR_cells <- "Others"
sc_obj$TCR_cells[rownames(sc_obj@meta.data) %in% TCRseq_barcode] <- "TCR Cells"
p4 = DimPlot(sc_obj, reduction = "umap", group.by = "TCR_cells", cols = c("grey50","blue"),pt.size = 0.2) + ggtitle("")
ggsave(filename = paste(output_dir, "P4.TCR.umap.png", sep="/"), plot = p4, width=9, height=7, path = "." )

# or 

# 提取UMAP坐标:
umap_data = sc_obj@reductions$umap@cell.embeddings %>% 
  as.data.frame() %>% rownames_to_column('barcode') %>%  
  tibble::column_to_rownames("barcode")
sc_obj = Seurat::AddMetaData(object=sc_obj, metadata=umap_data)

library(ggplot2)
# p4 = ggplot(data=metadata, aes(x=UMAP_1, y=UMAP_2, color=frequency,))+geom_point(size=log(metadata$frequency+1))+
p4 = ggplot(data=sc_obj@meta.data, aes(x=UMAP_1, y=UMAP_2, color=frequency,))+ geom_point(size=1)+
  scale_color_gradient(low="blue", high="red") +
  theme_classic()+
  theme(legend.title = element_blank(),
        panel.background = element_rect(fill = "transparent",colour = NA),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        plot.background = element_rect(fill = "transparent",colour = NA)
  )
p4
ggsave(filename =paste(output_dir, "p4.TCR.transparent.UMAP.png",sep = "/"), p4, device="png", width=1200, height = 900, units = "px",dpi = 72, bg="transparent")



########################## Find doublets
## doubletfinder
## install.packages('remotes')
## remotes::install_github('chris-mcginnis-ucsf/DoubletFinder')
library(DoubletFinder)
# Do not apply DoubletFinder to aggregated scRNA-seq data representing multiple distinct samples
## -------------- pK Identification (no ground-truth) --------------------------------------
find_doublets <- function(sc_obj){
  sweep.res.list <- DoubletFinder::paramSweep_v3(sc_obj, PCs = 1:10, sct = TRUE)
  sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)
  bcmvn <- find.pK(sweep.stats)
  mpK <- as.numeric(as.vector(bcmvn$pK[which.max(bcmvn$BCmetric)]))

  ## ----------------- DoubletFinder:Homotypic Doublet Proportion Estimate ------------------------
  annotations <- sc_obj@meta.data$seurat_clusters
  homotypic.prop <- modelHomotypic(annotations) 
  nExp_poi <- round(0.075*nrow(sc_obj@meta.data)) ## Assuming 7.5% doublet formation rate
  nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
  
  ## ------------DoubletFinder:Run DoubletFinder with varying classification stringencies ---------
  sc_obj <- doubletFinder_v3(sc_obj, PCs = 1:10, pN = 0.25, pK = mpK, 
                                    nExp = nExp_poi, reuse.pANN = FALSE, sct = TRUE)
  sc_obj <- doubletFinder_v3(sc_obj, PCs = 1:10, pN = 0.25, pK = mpK, 
                                    nExp = nExp_poi.adj, reuse.pANN = paste("pANN_0.25",mpK,nExp_poi,sep = "_"), sct = TRUE)
  
  ## ---------------------Doublets by DoubletFinder --------------------------------------------
  sc_obj$DF = sc_obj@meta.data[,match(paste("DF.classifications_0.25",mpK,nExp_poi,sep = "_"),colnames(sc_obj@meta.data))]
  doublets_rate = round(length(which(sc_obj$DF=="Doublet"))/ nrow(sc_obj@meta.data)*100,2)
  print(paste("doublets_rate : ", doublets_rate,"%", sep=""))
  sc_obj$DF[which((sc_obj$DF == "Doublet") & (sc_obj@meta.data[,match(paste("DF.classifications_0.25",mpK,nExp_poi.adj,sep = "_"),colnames(sc_obj@meta.data))] == "Singlet"))] <- "Doublet_lo"
  sc_obj$DF[which(sc_obj$DF == "Doublet")] <- "Doublet_hi"
  
  # write.table(sc_obj$DF, file = paste(output_dir, "P4.doublets.table.csv", sep="/"), sep = ",", row.names = FALSE, quote = FALSE)
  p5 = DimPlot(sc_obj, reduction = "umap", group.by = "DF", label.size = 5, label = TRUE, pt.size = 1) + 
  	ggtitle(paste("Doublets rate : ",doublets_rate,"%", sep=""))
  ggsave(filename = paste(output_dir, "P5.doublets.png", sep="/"), plot = p4, width=9, height=7, path = "./" )
  return(sc_obj)
}
sc_obj = find_doublets(sc_obj)
#########################

# make sure to choice these steps
# sc_obj = subset(sc_obj, subset=DF=='Singlet')
# sc_obj # 33538 features across 3283 samples within 1 assay 

# Using singR annotate cluster
## install.packages("BiocManager")
## BiocManager::install("SingleR")
## or 
## install.packages("devtools")
## devtools::install_github('dviraran/SingleR')

library(SingleR)
## BiocManager::install("celldex")
library(celldex)

sc_obj_SingleR = GetAssayData(sc_obj, slot="data")
sc_obj_clusters = sc_obj@meta.data$seurat_clusters

hpca = celldex::HumanPrimaryCellAtlasData()
blueprint_encode = celldex::BlueprintEncodeData()
dice = celldex::DatabaseImmuneCellExpressionData()
mona = celldex::MonacoImmuneData()
novershter = celldex::NovershternHematopoieticData()

pred.hpca = SingleR(test = sc_obj_SingleR, ref = hpca, labels = hpca$label.main,clusters = sc_obj_clusters )
pred.blue = SingleR(test = sc_obj_SingleR, ref = blueprint_encode, labels = blueprint_encode$label.main,clusters = sc_obj_clusters )
pred.dice = SingleR(test = sc_obj_SingleR, ref = dice, labels = dice$label.main,clusters = sc_obj_clusters)
pred.mona = SingleR(test = sc_obj_SingleR, ref = mona, labels = mona$label.main,clusters = sc_obj_clusters)
pred.novershter = SingleR(test = sc_obj_SingleR, ref = novershter, labels = novershter$label.main,clusters = sc_obj_clusters)

sc_obj_cellType = data.frame(ClusterID=levels(sc_obj@meta.data$seurat_clusters),
                                    hpca = pred.hpca$labels,
                                    blue  = pred.blue$labels,
                                    Dice = pred.dice$labels,
                                    mona = pred.mona$labels,
                                    novershter = pred.novershter$labels
)
write.table(sc_obj_cellType, file = paste(output_dir, "P6.SingleR_cell_type.csv",sep="/"),sep = ",", row.names = FALSE, quote=FALSE)
sc_obj@meta.data$singleR.hpca = sc_obj_cellType[match(sc_obj_clusters, sc_obj_cellType$ClusterID),'hpca']
p6_0 = DimPlot(sc_obj, reduction = "umap", label = T, group.by = 'singleR.hpca',)+ggtitle("HumanPrimaryCellAtlasData")
sc_obj@meta.data$singleR.blue = sc_obj_cellType[match(sc_obj_clusters, sc_obj_cellType$ClusterID),'blue']
p6_1 = DimPlot(sc_obj, reduction = "umap", label = T, group.by = 'singleR.blue',)+ggtitle("BlueprintEncodeData")
sc_obj@meta.data$singleR.dice = sc_obj_cellType[match(sc_obj_clusters, sc_obj_cellType$ClusterID),'Dice']
p6_2 = DimPlot(sc_obj, reduction = "umap", label = T, group.by = 'singleR.dice',)+ggtitle("DatabaseImmuneCellExpressionData")
sc_obj@meta.data$singleR.mona = sc_obj_cellType[match(sc_obj_clusters, sc_obj_cellType$ClusterID),'mona']
p6_3 = DimPlot(sc_obj, reduction = "umap", label = T, group.by = 'singleR.mona') + ggtitle("MonacoImmuneData")
sc_obj@meta.data$singleR.nove = sc_obj_cellType[match(sc_obj_clusters, sc_obj_cellType$ClusterID),'novershter']
p6_4 = DimPlot(sc_obj, reduction = "umap", label = T, group.by = 'singleR.nove',)+ ggtitle("NovershternHematopoieticData")
ggsave(filename=paste(output_dir,"P6_0.cluster.annotation.with.hpca.png",sep="/"),width=9,height=7,plot = p6_0, path = ".")
ggsave(filename=paste(output_dir,"P6_1.cluster.annotation.with.blue.png",sep="/"),width=9,height=7,plot = p6_1, path = ".")
ggsave(filename=paste(output_dir,"P6_2.cluster.annotation.with.dice.png",sep="/"),width=9,height=7,plot = p6_2, path = ".")
ggsave(filename=paste(output_dir,"P6_3.cluster.annotation.with.mona.png",sep="/"),width=9,height=7,plot = p6_3, path = ".")
ggsave(filename=paste(output_dir,"P6_4.cluster.annotation.with.nove.png",sep="/"),width=9,height=7,plot = p6_4, path = ".")

p6 = cowplot::plot_grid(p6_0,p6_1,p6_2,p6_3,p6_4, ncol = 2)
p6
ggsave(filename=paste(output_dir,"P6_cluster.annotation.with.5.databases.png",sep="/"),plot=p6,width=9,height=7, path = "./" )
## save.image(file=paste(output_dir, "sc_obj.RData",sep="/"), version = NULL, ascii = FALSE, safe = TRUE)

########################################
########################################
## 以下是个性化部分
# sc_obj@meta.data$singleR.target = sc_obj_cellType[match(sc_obj_clusters, sc_obj_cellType$ClusterID),'mona']
# sc_obj$singleR.target[which(sc_obj$clonotype_id == "clonotype597" )] <- "clon597is215"
# sc_obj$singleR.target[which(sc_obj$clonotype_id == "clonotype34" )] <- "clon34is1078"

# ump_data_for_ggplot  = Seurat::Embeddings(sc_obj, reduction = "umap")  %>%
#   as.data.frame() %>% rownames_to_column(var = 'barcode') %>%
#   left_join(data.frame("barcode"=rownames(sc_obj@meta.data),sc_obj@meta.data), by = 'barcode')
# ump_data_for_ggplot
# write.table(ump_data_for_ggplot, file=paste(output_dir, 'test.data.mass.csv', sep="/"), sep=",", row.names=FALSE)

# highlight_data <- ump_data_for_ggplot %>% filter(grepl("clon", singleR.target))
# highlight_data
# cells.highlight = rownames(data.frame("barcode"=rownames(sc_obj@meta.data),sc_obj@meta.data) %>% filter(grepl("clon", singleR.target)))
# cells.highlight
# p5_t = DimPlot(sc_obj, reduction = "umap", label = T, group.by = 'singleR.target',)+ ggtitle("MonacoImmuneData")
# ggsave(filename=paste(output_dir,"P5_t.cluster.annotation.with.targ.png",sep="/"),plot = p5_t,width=9,height=7, path = ".")
# p5_t

# or 
# 
# library(ggforce)
# install.packages("ggforce")
# library(tidyverse)
# library(magrittr)
# library(RColorBrewer)
# install.packages("ggrepel")
# library(ggrepel)
# brewer.pal(7,"Set2")[2]
# brewer.pal.info
# display.brewer.pal(11,"Paired")
# 
# tissue_colors <- c(`CD4+ T cells`= "red",
#                    `CD8+ T cells` = "green",
#                    `clonotype1004` = "#e491c1",
#                    `clonotype1012` = "#45a132",
#                    `clonotype1382` = "#45a132",
#                    `clonotype202` = "#6a6c2c",
#                    `clonotype245` = "#d8ab6a",
#                    `clonotype42` = "#e6c241",
#                    `clonotype215` = "#c44c90",
#                    `clonotype1078` = "#a9e082",
#                    `Dendritic cells` = "#f5899e",
#                    `T cells` = "blue"
# )
# cell_labels = c("CD4+ T cells", "CD8+ T cells", "clonotype1004", "clonotype1012", "clonotype1078", "clonotype1382", "clonotype202", "clonotype215", "clonotype42", "Dendritic cells", "T cells")
# 
# ump_data_for_ggplot  = Seurat::Embeddings(sc_obj, reduction = "umap")  %>% 
#   as.data.frame() %>% rownames_to_column(var = 'barcode') %>% 
#   left_join(data.frame("barcode"=rownames(sc_obj@meta.data),sc_obj@meta.data), by = 'barcode')
# ump_data_for_ggplot
# 
# highlight_data <- ump_data_for_ggplot %>% filter(grepl("clon", singleR.target))
# highlight_data
# common_themes = theme(plot.title = element_text(size = 17, hjust = 0.5),
#                       panel.grid.major = element_line(colour = NA), # 去掉网格线
#                       panel.background = element_blank(), # 去掉背景
#                       panel.grid.major.x = element_line(),  # 横向网格线linetype=2,xuxian
#                       panel.grid.major.y = element_line(), # 纵向网格线
#                       panel.border = element_rect(fill=NA), # 边框
#                       axis.line = element_line(), # 坐标轴线
#                       axis.title = element_text(size = 17),
#                       legend.title = element_blank(),
#                       #legend.position = c(0.95,0.95), # bottom, right, 'left'
#                       legend.background = element_rect(colour = NA, fill = NA),
# )
# p5_t2 = ggplot(ump_data_for_ggplot, aes(x=UMAP_1, y=UMAP_2, color=singleR.target, label=singleR.target))+ geom_point(size=0.5, alpha=1)+
#   scale_color_manual(values = tissue_colors) + geom_point(data = highlight_data, aes(x=UMAP_1, y=UMAP_2),size=2,) + 
#   geom_text_repel(data = highlight_data) +
#   common_themes
# 
# ggsave(filename=paste(output_dir,"P5_t2.cluster.annotation.with.targ.png",sep="/"),plot = p5_t2,width=9,height=7, path = ".")
# p5_t2
# 
# 
# ########################################
# ########################################
# ## sc_obj@meta.data # sc_obj[[]]
# ## write.table(sc_obj@meta.data, file=paste(output_dir, 'Total.meta.data.mass.csv', sep="/"), sep=",", row.names=TRUE, col.names = NA)
# write.table(data.frame("barcode"=rownames(sc_obj@meta.data),sc_obj@meta.data), file=paste(output_dir, 'Total.meta.data.mass.csv', sep="/"), sep=",", row.names=FALSE)

