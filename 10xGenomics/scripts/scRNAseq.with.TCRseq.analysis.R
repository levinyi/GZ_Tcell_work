# install.packages("Seurat")
# install.packages('tidyr')
# install.packages("tidyverse")
library(Seurat)
library(sctransform)
library(tidyverse, quietly=T)
library(cowplot)
args = commandArgs(T)

## setting for linux
# scRNAseq_path = args[1]
# tcr_path = args[2]
# project_name = args[3]

setwd("/cygene2/pipeline/10X/data/CR511_Post_infusion")
scRNAseq_path = "/cygene2/pipeline/10X/data/CR511_Post_infusion/CR511P02T-G507E6_GEX"
tcr_path = "/cygene2/pipeline/10X/data/CR511_Post_infusion/CR511P02T-G507E6_VDJ"

project_name = "CR511P02T-G507E6"

# create a output directory in current path.
output_dir = paste0(project_name,"-GEX-VDJ-analysis")
dir.create(output_dir)


# read scRNAseq folder:
data = Seurat::Read10X(data.dir = paste(scRNAseq_path,'outs/filtered_feature_bc_matrix', sep = "/"))

srt_obj <- Seurat::CreateSeuratObject(data, project = project_name)
srt_obj

load_vdj <- function(tcr_path){
  tcr_cells <- read.csv(paste(tcr_path, "outs/filtered_contig_annotations.csv", sep = "/"))
  # remove contigs that are not productive
  has_type <- grepl(pattern = "clonotype", x = tcr_cells$raw_clonotype_id)
  is_cell <- as.logical(tcr_cells$is_cell)
  confident <- as.logical(tcr_cells$high_confidence)
  is_fl <- as.logical(tcr_cells$full_length)
  is_productive <- as.logical(tcr_cells$productive)
  filter_tcrs <- has_type & is_cell & confident & is_fl & is_productive
  tcr_cells <- tcr_cells[filter_tcrs, ]
  
  total_counts <- length(unique(tcr_cells$barcode))
  tcr_cells[["total_counts"]] <- total_counts
  
  cdr3s_aa <- sapply(X = tcr_cells$raw_clonotype_id, FUN = generate_cdr3s, clone_info = tcr_cells)
  tcr_cells[["cdr3s_aa"]] <- cdr3s_aa
  
  clonalities <- sapply(X = tcr_cells$raw_clonotype_id,FUN = calc_clonality,clone_info = tcr_cells)
  tcr_cells[["TCR_clonality"]] <- clonalities
  
  v_genes <- sapply(X=tcr_cells$barcode, FUN = generate_genes, clone_info = tcr_cells, gene_type="v_gene")
  tcr_cells[["v_genes"]] <- v_genes
  j_genes <- sapply(X=tcr_cells$barcode, FUN = generate_genes, clone_info = tcr_cells, gene_type="j_gene")
  tcr_cells[["j_genes"]] <- j_genes
  
  
  clone_info <- as.data.frame(tcr_cells)
  return(clone_info)
}
generate_genes <- function(barcode, clone_info, gene_type) {
  rows <- which(clone_info$barcode == barcode)
  v_genes <- clone_info[rows, gene_type]
  return(v_genes)
}
generate_cdr3s <- function(clone, clone_info) {
  rows <- which(clone_info$raw_clonotype_id == clone)
  chains <- clone_info[rows, "chain"]
  cdr3s <- clone_info[rows, "cdr3"]
  cdr3_aa <- paste(chains, cdr3s, sep = ":") %>%
    unique %>% sort %>% paste(collapse = ";")
  return(cdr3_aa)
}

calc_clonality <- function(clone, clone_info) {
  rows <- which(clone_info$raw_clonotype_id == clone)
  clonality <- clone_info[rows, "barcode"] %>% unique %>% length
  return(clonality)
}

clone_info <- load_vdj(tcr_path = tcr_path)
head(clone_info)
tcr_data <- clone_info[!duplicated(clone_info$barcode),] %>% 
  dplyr::select(barcode, v_genes,j_genes,cdr3s_aa,TCR_clonality,raw_clonotype_id) 


tcr_data <- remove_rownames(tcr_data) %>% tibble::column_to_rownames("barcode")
head(tcr_data)
# read TCRseq file
srt_obj <- Seurat::AddMetaData(srt_obj, tcr_data)
head(srt_obj)

################
srt_obj[["percent.mt"]] <- PercentageFeatureSet(object = srt_obj, pattern = "^MT-")
srt_obj[["pct_ribo"]]   <- PercentageFeatureSet(object = srt_obj, pattern = "^RP[SL]")

VlnPlot(srt_obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "pct_ribo"), ncol = 2)
# ggsave(filename=paste(output_dir,"P1.VlnPlot.feature.png",sep="/"), plot=p1, width=12,height=6, path=".")
# FeatureScatter is typically used to visualize feature-feature relationships, but can be used
# for anything calculated by the object, i.e. columns in object metadata, PC scores etc.
head(srt_obj@meta.data)


# plot1 <- FeatureScatter(srt_obj, feature1 = "nCount_RNA", feature2 = "percent.mt")
# plot2 <- FeatureScatter(srt_obj, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
# p2 = CombinePlots(plots = list(plot1, plot2))
# ggsave(filename=paste(output_dir, "P2.combinePlots.FeatureScatter.png", sep="/"), plot=p2,width=12,height=6, path=".")
ggplot(data = srt_obj@meta.data, aes(x=percent.mt)) + geom_density(position = "stack")

srt_obj <- subset(srt_obj, subset = nFeature_RNA > 20 & nFeature_RNA < 3000 & percent.mt < 25)
# srt_obj # 33538 features across 3549 samples within 1 assay 
# standard seurat workflow
srt_obj <- SCTransform(srt_obj, method = "glmGamPoi", vars.to.regress = "percent.mt", verbose = FALSE) %>%
  RunPCA(npcs=30) %>%
  FindNeighbors(dims = 1:10) %>%
  RunUMAP(dims = 1:10) %>%
  FindClusters(resolution=0.5) 

p3 <- DimPlot(srt_obj, reduction = "umap")
p3
ggsave(filename = paste(output_dir, "P3.raw.umap.png", sep="/"), plot = p3, width=9, height=7, path = "." )


FeaturePlot(srt_obj, features = c("CD4","CD8A","CD8B"))
###########################################################
### for TCR data
TCRseq_barcode = rownames(srt_obj@meta.data[which(srt_obj@meta.data$clonotype_id != "NA"),])
length(TCRseq_barcode)

# DimPlot(srt_obj, reduction = "umap", cells.highlight = TCRseq_barcode, pt.size = 0.3)+
#   scale_color_manual(labels = c("Others","TCR cells"),values = c("grey50","blue")) #  + labs(color = "legend title")
# or
srt_obj$TCR_cells <- "Others"
srt_obj$TCR_cells[rownames(srt_obj@meta.data) %in% TCRseq_barcode] <- "TCR Cells"
p4 = DimPlot(srt_obj, reduction = "umap", group.by = "TCR_cells", cols = c("grey50","blue"),pt.size = 0.2) + ggtitle("")
ggsave(filename = paste(output_dir, "P4.TCR.umap.png", sep="/"), plot = p4, width=9, height=7, path = "." )

# or 

# 提取UMAP坐标:
umap_data = srt_obj@reductions$umap@cell.embeddings %>% 
  as.data.frame() %>% rownames_to_column('barcode') %>%  
  tibble::column_to_rownames("barcode")
srt_obj = Seurat::AddMetaData(object=srt_obj, metadata=umap_data)

library(ggplot2)
# p4 = ggplot(data=metadata, aes(x=UMAP_1, y=UMAP_2, color=frequency,))+geom_point(size=log(metadata$frequency+1))+
p4 = ggplot(data=srt_obj@meta.data, aes(x=UMAP_1, y=UMAP_2, color=frequency,))+ geom_point(size=1)+
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
find_doublets <- function(srt_obj){
  sweep.res.list <- DoubletFinder::paramSweep_v3(srt_obj, PCs = 1:10, sct = TRUE)
  sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)
  bcmvn <- find.pK(sweep.stats)
  mpK <- as.numeric(as.vector(bcmvn$pK[which.max(bcmvn$BCmetric)]))

  ## ----------------- DoubletFinder:Homotypic Doublet Proportion Estimate ------------------------
  annotations <- srt_obj@meta.data$seurat_clusters
  homotypic.prop <- modelHomotypic(annotations) 
  nExp_poi <- round(0.075*nrow(srt_obj@meta.data)) ## Assuming 7.5% doublet formation rate
  nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
  
  ## ------------DoubletFinder:Run DoubletFinder with varying classification stringencies ---------
  srt_obj <- doubletFinder_v3(srt_obj, PCs = 1:10, pN = 0.25, pK = mpK, 
                                    nExp = nExp_poi, reuse.pANN = FALSE, sct = TRUE)
  srt_obj <- doubletFinder_v3(srt_obj, PCs = 1:10, pN = 0.25, pK = mpK, 
                                    nExp = nExp_poi.adj, reuse.pANN = paste("pANN_0.25",mpK,nExp_poi,sep = "_"), sct = TRUE)
  
  ## ---------------------Doublets by DoubletFinder --------------------------------------------
  srt_obj$DF = srt_obj@meta.data[,match(paste("DF.classifications_0.25",mpK,nExp_poi,sep = "_"),colnames(srt_obj@meta.data))]
  doublets_rate = round(length(which(srt_obj$DF=="Doublet"))/ nrow(srt_obj@meta.data)*100,2)
  print(paste("doublets_rate : ", doublets_rate,"%", sep=""))
  srt_obj$DF[which((srt_obj$DF == "Doublet") & (srt_obj@meta.data[,match(paste("DF.classifications_0.25",mpK,nExp_poi.adj,sep = "_"),colnames(srt_obj@meta.data))] == "Singlet"))] <- "Doublet_lo"
  srt_obj$DF[which(srt_obj$DF == "Doublet")] <- "Doublet_hi"
  
  # write.table(srt_obj$DF, file = paste(output_dir, "P4.doublets.table.csv", sep="/"), sep = ",", row.names = FALSE, quote = FALSE)
  p5 = DimPlot(srt_obj, reduction = "umap", group.by = "DF", label.size = 5, label = TRUE, pt.size = 1) + 
  	ggtitle(paste("Doublets rate : ",doublets_rate,"%", sep=""))
  ggsave(filename = paste(output_dir, "P5.doublets.png", sep="/"), plot = p4, width=9, height=7, path = "./" )
  return(srt_obj)
}
srt_obj = find_doublets(srt_obj)
#########################

# make sure to choice these steps
# srt_obj = subset(srt_obj, subset=DF=='Singlet')
# srt_obj # 33538 features across 3283 samples within 1 assay 

# Using singR annotate cluster
## install.packages("BiocManager")
## BiocManager::install("SingleR")
## or 
## install.packages("devtools")
## devtools::install_github('dviraran/SingleR')

library(SingleR)
## BiocManager::install("celldex")
library(celldex)

srt_obj_SingleR = GetAssayData(srt_obj, slot="data")
srt_obj_clusters = srt_obj@meta.data$seurat_clusters

hpca = celldex::HumanPrimaryCellAtlasData()
blueprint_encode = celldex::BlueprintEncodeData()
dice = celldex::DatabaseImmuneCellExpressionData()
mona = celldex::MonacoImmuneData()
novershter = celldex::NovershternHematopoieticData()

pred.hpca = SingleR(test = srt_obj_SingleR, ref = hpca, labels = hpca$label.main,clusters = srt_obj_clusters )
pred.blue = SingleR(test = srt_obj_SingleR, ref = blueprint_encode, labels = blueprint_encode$label.main,clusters = srt_obj_clusters )
pred.dice = SingleR(test = srt_obj_SingleR, ref = dice, labels = dice$label.main,clusters = srt_obj_clusters)
pred.mona = SingleR(test = srt_obj_SingleR, ref = mona, labels = mona$label.main,clusters = srt_obj_clusters)
pred.novershter = SingleR(test = srt_obj_SingleR, ref = novershter, labels = novershter$label.main,clusters = srt_obj_clusters)

srt_obj_cellType = data.frame(ClusterID=levels(srt_obj@meta.data$seurat_clusters),
                                    hpca = pred.hpca$labels,
                                    blue  = pred.blue$labels,
                                    Dice = pred.dice$labels,
                                    mona = pred.mona$labels,
                                    novershter = pred.novershter$labels
)
write.table(srt_obj_cellType, file = paste(output_dir, "P6.SingleR_cell_type.csv",sep="/"),sep = ",", row.names = FALSE, quote=FALSE)
srt_obj@meta.data$singleR.hpca = srt_obj_cellType[match(srt_obj_clusters, srt_obj_cellType$ClusterID),'hpca']
p6_0 = DimPlot(srt_obj, reduction = "umap", label = T, group.by = 'singleR.hpca',)+ggtitle("HumanPrimaryCellAtlasData")
srt_obj@meta.data$singleR.blue = srt_obj_cellType[match(srt_obj_clusters, srt_obj_cellType$ClusterID),'blue']
p6_1 = DimPlot(srt_obj, reduction = "umap", label = T, group.by = 'singleR.blue',)+ggtitle("BlueprintEncodeData")
srt_obj@meta.data$singleR.dice = srt_obj_cellType[match(srt_obj_clusters, srt_obj_cellType$ClusterID),'Dice']
p6_2 = DimPlot(srt_obj, reduction = "umap", label = T, group.by = 'singleR.dice',)+ggtitle("DatabaseImmuneCellExpressionData")
srt_obj@meta.data$singleR.mona = srt_obj_cellType[match(srt_obj_clusters, srt_obj_cellType$ClusterID),'mona']
p6_3 = DimPlot(srt_obj, reduction = "umap", label = T, group.by = 'singleR.mona') + ggtitle("MonacoImmuneData")
srt_obj@meta.data$singleR.nove = srt_obj_cellType[match(srt_obj_clusters, srt_obj_cellType$ClusterID),'novershter']
p6_4 = DimPlot(srt_obj, reduction = "umap", label = T, group.by = 'singleR.nove',)+ ggtitle("NovershternHematopoieticData")
ggsave(filename=paste(output_dir,"P6_0.cluster.annotation.with.hpca.png",sep="/"),width=9,height=7,plot = p6_0, path = ".")
ggsave(filename=paste(output_dir,"P6_1.cluster.annotation.with.blue.png",sep="/"),width=9,height=7,plot = p6_1, path = ".")
ggsave(filename=paste(output_dir,"P6_2.cluster.annotation.with.dice.png",sep="/"),width=9,height=7,plot = p6_2, path = ".")
ggsave(filename=paste(output_dir,"P6_3.cluster.annotation.with.mona.png",sep="/"),width=9,height=7,plot = p6_3, path = ".")
ggsave(filename=paste(output_dir,"P6_4.cluster.annotation.with.nove.png",sep="/"),width=9,height=7,plot = p6_4, path = ".")

p6 = cowplot::plot_grid(p6_0,p6_1,p6_2,p6_3,p6_4, ncol = 2)
p6
ggsave(filename=paste(output_dir,"P6_cluster.annotation.with.5.databases.png",sep="/"),plot=p6,width=9,height=7, path = "./" )
## save.image(file=paste(output_dir, "srt_obj.RData",sep="/"), version = NULL, ascii = FALSE, safe = TRUE)

########################################
########################################
## 以下是个性化部分
# srt_obj@meta.data$singleR.target = srt_obj_cellType[match(srt_obj_clusters, srt_obj_cellType$ClusterID),'mona']
# srt_obj$singleR.target[which(srt_obj$clonotype_id == "clonotype597" )] <- "clon597is215"
# srt_obj$singleR.target[which(srt_obj$clonotype_id == "clonotype34" )] <- "clon34is1078"

# ump_data_for_ggplot  = Seurat::Embeddings(srt_obj, reduction = "umap")  %>%
#   as.data.frame() %>% rownames_to_column(var = 'barcode') %>%
#   left_join(data.frame("barcode"=rownames(srt_obj@meta.data),srt_obj@meta.data), by = 'barcode')
# ump_data_for_ggplot
# write.table(ump_data_for_ggplot, file=paste(output_dir, 'test.data.mass.csv', sep="/"), sep=",", row.names=FALSE)

# highlight_data <- ump_data_for_ggplot %>% filter(grepl("clon", singleR.target))
# highlight_data
# cells.highlight = rownames(data.frame("barcode"=rownames(srt_obj@meta.data),srt_obj@meta.data) %>% filter(grepl("clon", singleR.target)))
# cells.highlight
# p5_t = DimPlot(srt_obj, reduction = "umap", label = T, group.by = 'singleR.target',)+ ggtitle("MonacoImmuneData")
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
# ump_data_for_ggplot  = Seurat::Embeddings(srt_obj, reduction = "umap")  %>% 
#   as.data.frame() %>% rownames_to_column(var = 'barcode') %>% 
#   left_join(data.frame("barcode"=rownames(srt_obj@meta.data),srt_obj@meta.data), by = 'barcode')
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
# ## srt_obj@meta.data # srt_obj[[]]
# ## write.table(srt_obj@meta.data, file=paste(output_dir, 'Total.meta.data.mass.csv', sep="/"), sep=",", row.names=TRUE, col.names = NA)
# write.table(data.frame("barcode"=rownames(srt_obj@meta.data),srt_obj@meta.data), file=paste(output_dir, 'Total.meta.data.mass.csv', sep="/"), sep=",", row.names=FALSE)

