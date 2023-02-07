#' This is easy function that can annotate scRNAseq object by just pass a object
#' to generate plot and table for annotation using SingleR package.
#' only for human species.


library(SingleR)
# BiocManager::install("celldex")
library(celldex)
args = commandArgs(T)


srt_obj <- readRDS(args[1])
srt_obj_SingleR = GetAssayData(srt_obj, slot="data")
srt_obj_clusters = srt_obj@meta.data$seurat_clusters


##### load from local.
# hpca = celldex::HumanPrimaryCellAtlasData()
# blueprint_encode = celldex::BlueprintEncodeData()
# dice = celldex::DatabaseImmuneCellExpressionData()
# mona = celldex::MonacoImmuneData()
# nove = celldex::novenHematopoieticData()
hpca = readRDS("xxx.rds")
blue = readRDS("xxx.rds")
dice = readRDS("xxx.rds")
mona = readRDS("xxx.rds")
nove = readRDS("xxx.rds")



# choose label.main or label.fine
pred.hpca = SingleR(test = srt_obj_SingleR, ref = hpca, labels = hpca$label.fine, clusters = srt_obj_clusters)
pred.blue = SingleR(test = srt_obj_SingleR, ref = blue, labels = blue$label.fine, clusters = srt_obj_clusters)
pred.dice = SingleR(test = srt_obj_SingleR, ref = dice, labels = dice$label.fine, clusters = srt_obj_clusters)
pred.mona = SingleR(test = srt_obj_SingleR, ref = mona, labels = mona$label.fine, clusters = srt_obj_clusters)
pred.nove = SingleR(test = srt_obj_SingleR, ref = nove, labels = nove$label.fine, clusters = srt_obj_clusters)

srt_obj_cellType = data.frame(ClusterID=levels(sc_seurat_obj@meta.data$seurat_clusters),
                                    hpca = pred.hpca$labels,
                                    blue = pred.blue$labels,
                                    Dice = pred.dice$labels,
                                    mona = pred.mona$labels,
                                    nove = pred.nove$labels
)
write.table(srt_obj_cellType, file = paste(output_dir,"P5.SingleR_cell_type.csv",sep="/"),sep = ",", row.names = FALSE, quote=FALSE)

sc_seurat_obj@meta.data$singleR.hpca = srt_obj_cellType[match(srt_obj_clusters, srt_obj_cellType$ClusterID),'hpca']
sc_seurat_obj@meta.data$singleR.blue = srt_obj_cellType[match(srt_obj_clusters, srt_obj_cellType$ClusterID),'blue']
sc_seurat_obj@meta.data$singleR.dice = srt_obj_cellType[match(srt_obj_clusters, srt_obj_cellType$ClusterID),'Dice']
sc_seurat_obj@meta.data$singleR.mona = srt_obj_cellType[match(srt_obj_clusters, srt_obj_cellType$ClusterID),'mona']
sc_seurat_obj@meta.data$singleR.nove = srt_obj_cellType[match(srt_obj_clusters, srt_obj_cellType$ClusterID),'nove']

p5_0 = DimPlot(sc_seurat_obj, reduction = "umap", label = T, group.by = 'singleR.hpca',)+ ggtitle("HumanPrimaryCellAtlasData")
p5_1 = DimPlot(sc_seurat_obj, reduction = "umap", label = T, group.by = 'singleR.blue',)+ ggtitle("BlueprintEncodeData")
p5_2 = DimPlot(sc_seurat_obj, reduction = "umap", label = T, group.by = 'singleR.dice',)+ ggtitle("DatabaseImmuneCellExpressionData")
p5_3 = DimPlot(sc_seurat_obj, reduction = "umap", label = T, group.by = 'singleR.mona') + ggtitle("MonacoImmuneData")
p5_4 = DimPlot(sc_seurat_obj, reduction = "umap", label = T, group.by = 'singleR.nove',)+ ggtitle("novenHematopoieticData")
ggsave(filename=paste(output_dir,"4Annotation/SingleR/P5_0.cluster.annotation.with.hpca.png",sep="/"),width=9,height=7,plot = p5_0, path = "./")
ggsave(filename=paste(output_dir,"4Annotation/SingleR/P5_1.cluster.annotation.with.blue.png",sep="/"),width=9,height=7,plot = p5_1, path = "./")
ggsave(filename=paste(output_dir,"4Annotation/SingleR/P5_2.cluster.annotation.with.dice.png",sep="/"),width=9,height=7,plot = p5_2, path = "./")
ggsave(filename=paste(output_dir,"4Annotation/SingleR/P5_3.cluster.annotation.with.mona.png",sep="/"),width=9,height=7,plot = p5_3, path = "./")
ggsave(filename=paste(output_dir,"4Annotation/SingleR/P5_4.cluster.annotation.with.nove.png",sep="/"),width=9,height=7,plot = p5_4, path = "./")
library(cowplot)
p5 = cowplot::plot_grid(p5_0,p5_1,p5_2,p5_3,p5_4, ncol = 2)
# p5 = cowplot::plot_grid(p5_0,p5_1,p5_3,p5_4, ncol = 2)
ggsave(filename=paste(output_dir,"4Annotation/SingleR/P5_cluster.annotation.with.5.databases.png",sep="/"), plot = p5, width=9, height=7, path = "./" )


