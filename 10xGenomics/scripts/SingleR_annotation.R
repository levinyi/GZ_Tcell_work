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
write.table(sc_seurat_obj_cellType, file = paste(output_dir,"P5.SingleR_cell_type.csv",sep="/"),sep = ",", row.names = FALSE, quote=FALSE)
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
ggsave(filename=paste(output_dir,"4Annotation/SingleR/P5_0.cluster.annotation.with.hpca.png",sep="/"),width=9,height=7,plot = p5_0, path = "./")
ggsave(filename=paste(output_dir,"4Annotation/SingleR/P5_1.cluster.annotation.with.blue.png",sep="/"),width=9,height=7,plot = p5_1, path = "./")
ggsave(filename=paste(output_dir,"4Annotation/SingleR/P5_2.cluster.annotation.with.dice.png",sep="/"),width=9,height=7,plot = p5_2, path = "./")
ggsave(filename=paste(output_dir,"4Annotation/SingleR/P5_3.cluster.annotation.with.mona.png",sep="/"),width=9,height=7,plot = p5_3, path = "./")
ggsave(filename=paste(output_dir,"4Annotation/SingleR/P5_4.cluster.annotation.with.nove.png",sep="/"),width=9,height=7,plot = p5_4, path = "./")
library(cowplot)
p5 = cowplot::plot_grid(p5_0,p5_1,p5_2,p5_3,p5_4, ncol = 2)
# p5 = cowplot::plot_grid(p5_0,p5_1,p5_3,p5_4, ncol = 2)
ggsave(filename=paste(output_dir,"4Annotation/SingleR/P5_cluster.annotation.with.5.databases.png",sep="/"), plot = p5, width=9, height=7, path = "./" )


