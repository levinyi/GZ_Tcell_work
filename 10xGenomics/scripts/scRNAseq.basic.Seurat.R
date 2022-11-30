library(Seurat)
library(sctransform)
library(glmGamPoi)
library(tidyverse, quietly=T)
args = commandArgs(T)

## setting for linux
scRNAseq_path = args[1]
# setwd("/cygene2/pipeline/10X/data/SA501001-TD/SA501001-G517E5_GEX")
# scRNAseq_path = "/cygene2/pipeline/10X/data/SA501001-TD/SA501001-G517E5_GEX"

output_dir = args[2]
# output_dir = "G517E5_GEX_basic_Seurat"

# read scRNAseq folder:
data = Seurat::Read10X(data.dir = paste(scRNAseq_path, 'outs/filtered_feature_bc_matrix', sep = "/"))

# output_RNA_data = data %>% as.matrix() %>% t() %>%  as.data.frame() %>%   rownames_to_column("barcode") %>% select(c("barcode","CD4","CD8A","CD8B","FOXP3"))
# write.table(output_RNA_data, file=paste(output_dir,'1QC/RNA.t.csv',sep="/"), sep = ",", row.names = FALSE)


srt_obj <- CreateSeuratObject(data) %>%
	PercentageFeatureSet(pattern = "^MT-", col.name = "percent.mt") %>%
	PercentageFeatureSet(pattern = "^RP[SL]", col.name = "pct_ribo")

p1 = VlnPlot(srt_obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "pct_ribo"), 
        log=TRUE, flip=T, ncol=2)
ggsave(filename=paste(output_dir,"P1.VlnPlot.feature.png",sep="/"), plot=p1, width=12, height=6)

srt_obj <- subset(srt_obj, subset = percent.mt <10 & pct_ribo >6)

srt_obj <-	SCTransform(srt_obj, method = "glmGamPoi", vars.to.regress = "percent.mt", verbose = FALSE) %>%
	RunPCA() %>%
	FindNeighbors(dims = 1:30) %>%
	RunUMAP(dims = 1:30)



# plot1 <- FeatureScatter(srt_obj, feature1 = "nCount_RNA", feature2 = "percent.mt")
# plot2 <- FeatureScatter(srt_obj, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
# p2 = CombinePlots(plots = list(plot1, plot2))
# ggsave(filename=paste(output_dir,"1QC/P2.combinePlots.FeatureScatter.png",sep="/"), plot=p2, width=12,height=6,path=".")


srt_obj <- FindClusters(srt_obj, resolution = c(0.4, 0.8, 1.2, 1.6, 2))
p3 <- DimPlot(srt_obj, reduction = "umap",label = T, group.by = "SCT_snn_res.0.4") +
    DimPlot(srt_obj, reduction = "umap",label = T, group.by = "SCT_snn_res.0.8") +
    DimPlot(srt_obj, reduction = "umap",label = T, group.by = "SCT_snn_res.1.2") +
    DimPlot(srt_obj, reduction = "umap",label = T, group.by = "SCT_snn_res.1.6") +
    DimPlot(srt_obj, reduction = "umap",label = T, group.by = "SCT_snn_res.2")
ggsave(filename = paste(output_dir,"GEX_UMAP_GexClusterByRes.png",sep="/"),plot = p3, width=18, height=12)


multi_feature_plot <- function(srt_obj, feature_list, ggplot_mods, 
                               point_size=0.5, order=FALSE,min_cut=NA,
                               max_cut=NA) {
    plot_list <- Seurat::FeaturePlot(srt_obj,features = feature_list, pt.size = point_size,combine = FALSE,order = order,min.cutoff = min_cut,max.cutoff = max_cut)
    plot_list <- lapply(X = plot_list, 
                        FUN = function(x) {
                            for (mod in ggplot_mods) {
                                x <- x + mod
                            }
                            return(x)
                        }
    )
    combined_plot <- Reduce(`+`, plot_list)
    return(combined_plot)
}


rootpath_jet_scaling <- function(jet_colors=RP_DS_JET,
                                 value_range=NULL, 
                                 null_color= ROOTPATH_COLOURS["black"],
                                 log=FALSE){
    transformation="identity"
    if (log) {
        transformation <- "log"
    }
    new_scale <- ggplot2::scale_color_gradientn(
        colours = jet_colors,limits=value_range, 
        na.value = null_color, trans=transformation
    )
    return(new_scale)
}
ROOTPATH_COLOURS <- list("black" = "#3E3D57", "violet" = "#FF00F1",
                         "orange" = "#FF996A", "yellow" = "#FFEC00",
                         "green" = "#00F300", "navy" = "#3565AA",
                         "blue" = "#0069FF", "purple" = "#7D2DFD",
                         "red" = "#FC5B68", "skyblue" = "#5FB3D3")
RP_DS_JET <- c(ROOTPATH_COLOURS["navy"], ROOTPATH_COLOURS["blue"],
               ROOTPATH_COLOURS["skyblue"], ROOTPATH_COLOURS["yellow"],
               ROOTPATH_COLOURS["orange"], ROOTPATH_COLOURS["red"])

yost_markers = c("CD3D","CD4","CD8A","CD8B","FOXP3","CCR7", "IL26", "CD200","EOMES","KLRD1","IFNG","HAVCR2","CXCL13","ENTPD1","PDCD1","ITGAE")
yost_markers <- yost_markers[yost_markers %in% rownames(srt_obj)]
yost_fp <- multi_feature_plot(srt_obj=srt_obj, feature_list = yost_markers, 
                              ggplot_mods=c(rootpath_jet_scaling())
)
ggsave(filename = paste(output_dir, "GEX_UMAP_YostFP.png",sep="/"),plot = yost_fp, width = 18, height = 12)