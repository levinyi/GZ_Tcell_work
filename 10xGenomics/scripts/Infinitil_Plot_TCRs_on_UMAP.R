#' example: Rscript  /cygene/work/00.test/pipeline/10xGenomics/scripts/Infinitil_Plot_TCRs_on_UMAP.R integrated_seurat.rds tcr_id.txt
#' \input tcr_id.txt
#' \input integrated_seurat.rds
#' \output two pdfs.

library(Seurat)
library(ggplot2)
library(magrittr)
library(tibble)

args = commandArgs(T)

seurat_rds = args[1] # "integrated_seurat.rds"
# seurat_rds = "/cygene2/pipeline/10X/data/CR511_Post_infusion/analysis/integrated_seurat.rds"
srt_obj = readRDS(seurat_rds)

tcr_file = args[2]
# tcr_file = "/cygene2/pipeline/10X/data/CR511_Post_infusion/CR511P02T_EntrySeq/endoTCR.bc.counts.txt"
df = read.table(tcr_file, header=T)

# check if there are duplicated values in the barcode column.
if (any(duplicated(df[,"barcode"]))){
    # if there are duplicated values, remove them and keep other columns
    df <- df[!duplicated(df[,"barcode"]),]
    rownames(df) <- NULL
}
tcr_table = df %>% tibble::column_to_rownames("barcode")
head(tcr_table)
# add metadata
srt_obj <- Seurat::AddMetaData(srt_obj, tcr_table)

# For p1, plot highlight cells by dimplot
# highlighted_barcode <- rownames(srt_obj@meta.data[srt_obj@meta.data$clonotype_id %in% tcr_id,])
highlighted_barcode <- rownames(tcr_table)
## 
p <- DimPlot(srt_obj, cells.highlight = highlighted_barcode, split.by = "orig.ident")
p <- p + scale_color_manual(
		values=c("grey50","red"),
		breaks=c("Unselected","Group_1"),
		labels=c("Others","Selected")
	)
p
ggsave("Plot.Selected.TCRs.on.UMAP.by.sample.pdf",plot=p, device="pdf",width=15, height=5)
ggsave("Plot.Selected.TCRs.on.UMAP.by.sample.png",plot=p, device="png",width=15, height=5)


# # for individual color and split sample
p_test <- DimPlot(srt_obj, group.by = "clonotype_ID", split.by = "orig.ident") + theme(plot.title = element_blank())
ggsave("Plot.Selected.TCRs.on.UMAP.by.sample.DimPlot.pdf",plot=p_test, device="pdf",width=15, height=5)
ggsave("Plot.Selected.TCRs.on.UMAP.by.sample.DimPlot.png",plot=p_test, device="png",width=15, height=5)


#########################################
# for individual color and split sample
data = FetchData(srt_obj, vars=c("UMAP_1", "UMAP_2", "orig.ident", "clonotype_ID"))

split_df <- split(data, data$orig.ident)
plot_list <- lapply(seq_along(split_df), function(i) {
    x = split_df[[i]]
    name= (names(split_df[i]))
    data_list = split(x, is.na(x$clonotype_ID))
    ggplot(data_list$`TRUE`, aes(x=UMAP_1, y=UMAP_2)) +
        geom_point(size=0.5,color="grey") +
        geom_point(data=data_list$`FALSE`, aes(x=UMAP_1, y=UMAP_2, color=clonotype_ID),size=1) +
        theme_classic() + theme(legend.title=element_blank(),legend.position = "none") +
        ylab("") +
        ggtitle(name) +
        theme(
            plot.title = element_text(hjust = 0.5),
            axis.ticks.y = element_blank(),
            axis.line.y = element_blank(),
            axis.text.y = element_blank())
})

# 修改最后一张图的图例
plot_list[[length(plot_list)]] <- plot_list[[length(plot_list)]] + 
    theme(legend.position = "right")
# 修改第一张图的图例
plot_list[[1]] <- plot_list[[1]] +theme(
    axis.ticks.y = element_line(),
    axis.line.y = element_line(),
    axis.text.y = element_text(),
    axis.title.y = element_text(),
) + ylab("UMP_2")
combined_plot <- Reduce(`+`, plot_list)
combined_plot

ggsave("Plot.Selected.TCRs.on.UMAP.sep.cor.by.sample.pdf", plot=combined_plot, device="pdf", width=15,height=5)
ggsave("Plot.Selected.TCRs.on.UMAP.sep.cor.by.sample.png", plot=combined_plot, device="png", width=15,height=5)


######################################################################################
# head(split_df)
# data$plot <- ifelse(is.na(data$clonotype_ID), "others", data$clonotype_ID)
# 
# # split selected data and Unselected data. give ggplot2 two layers.
# selected_data = data %>% dplyr::filter(plot != "others")
# unselect_data = data %>% dplyr::filter(plot == "others")
# 
# # ggplot2
# p2 <- ggplot(unselect_data, aes(x=UMAP_1,y=UMAP_2)) +
# 	geom_point(size=0.7, alpha=1, color="grey") +
# 	geom_point(data=selected_data, aes(x=UMAP_1,y=UMAP_2, color=clonotype_ID),size=1) +
# 	theme_classic() + theme(legend.title=element_blank())
# p2

# convert group1_clonotype80 to clonotype80: 
###short_name = c()
###for (each in tcr_id) {
###	if (!startsWith(each, "clonotype")){
###		print(paste("this clonotype is not starts with clonotype", each, sep=":"))
###		name = substring(each, 8)
###		short_name = c(short_name, name)
###	}
###}
###print(short_name)
###p2 <- p2 + scale_color_discrete(breaks=tcr_id, labels=short_name)
# ggsave("Plot.Selected.TCRs.on.UMAP.sep.cor.pdf", plot=p2, device="pdf", width=7.5,height=5)
# ggsave("Plot.Selected.TCRs.on.UMAP.sep.cor.png", plot=p2, device="png", width=7.5,height=5)
# 
# ##############################################
# # for counts plot
# src("Infinitil_Screening_analysis.R")
# plot_spots(srt_obj, features="counts", save_path=save_path)
# print("done")

#########################################################
#########################################################Data to Drag
FeaturePlot(srt_obj, features = "counts",split.by = "orig.ident")
RP_GRADIENT <- c(ROOTPATH_COLOURS["navy"],
                 ROOTPATH_COLOURS["yellow"],
                 ROOTPATH_COLOURS["violet"])
ggplot_mods <- c(rootpath_jet_scaling())
plot_spots <- function(srt_obj, features, save_path){
    plot_list <- Seurat::FeaturePlot(object = srt_obj,
                                features = features, order = TRUE,
                                split.by = "orig.ident",
                                combine = F)
    # plot_list <- lapply(X=plot_list, FUN = function(x){
    #     x <- x + rootpath_jet_scaling(RP_GRADIENT, log=TRUE)
    # })
    # plot = Reduce(`+`, plot_list)
    # return(plot)
    plot_list <- lapply(X = plot_list,FUN = function(x) {
        for (mod in ggplot_mods) {
            x <- x + mod
        }
        return(x)
    })
    plot <- Reduce(`+`, plot_list)
    return(plot)
}
plot = plot_spots(srt_obj, features="counts")
plot
ggsave("Plot.BarcodeCounts.on.UMAP.sep.cor.by.sample.pdf", plot=plot, device="pdf", width=15,height=5)
ggsave("Plot.BarcodeCounts.on.UMAP.sep.cor.by.sample.png", plot=plot, device="png", width=15,height=5)

##################################################################
### method 2 by ggplot2

data = FetchData(srt_obj, vars=c("UMAP_1","UMAP_2","counts","orig.ident"))
head(data)
split_df <- split(data, data$orig.ident)
# color
# size = counts
RP_DS_JET <- c(ROOTPATH_COLOURS["navy"], ROOTPATH_COLOURS["blue"],
               ROOTPATH_COLOURS["skyblue"], ROOTPATH_COLOURS["yellow"],
               ROOTPATH_COLOURS["orange"], ROOTPATH_COLOURS["red"])

ROOTPATH_COLOURS <- list("black" = "#3E3D57", "violet" = "#FF00F1",
                         "orange" = "#FF996A", "yellow" = "#FFEC00",
                         "green" = "#00F300", "navy" = "#3565AA",
                         "blue" = "#0069FF", "purple" = "#7D2DFD",
                         "red" = "#FC5B68", "skyblue" = "#5FB3D3")

rootpath_jet_scaling <- function(jet_colors = RP_DS_JET,
                                 value_range = NULL,
                                 null_color = ROOTPATH_COLOURS["black"],
                                 log = FALSE
) {
    transformation = "identity"
    if (log) {
        transformation <- "log"
    }
    new_scale <- ggplot2::scale_color_gradientn(colours = jet_colors,
                                                limits = value_range,
                                                na.value = null_color,
                                                trans = transformation
    )
    return(new_scale)
}


ggplot_mods <- c(rootpath_jet_scaling())
ggplot_mods <- c(rootpath_jet_scaling(value_range = c(0, 0.25)))
ggplot_mods

plot_list <- lapply(seq_along(split_df), function(i){
    x = split_df[[i]]
    name = names(split_df[i])
    # print(name)
    data_list <- split(x, is.na(x$counts))
    ggplot(data_list$`TRUE`,aes(UMAP_1,UMAP_2)) + geom_point(color="grey",size=0.5) +
        geom_point(data=data_list$`FALSE`, aes(UMAP_1,UMAP_2),size=data_list$`FALSE`$counts) +
        theme_classic() + theme(legend.title = element_blank(), legend.position = "none") +
        ylab("") + ggtitle(name)+ theme(
            plot.title = element_text(hjust = 0.5),
            axis.ticks.y = element_blank(),
            axis.line.y = element_blank(),
            axis.text.y = element_blank()
        ) + ggplot_mods
})
# add legend
plot_list[[length(plot_list)]] <- plot_list[[length(plot_list)]] + theme(legend.position = "right")
# 
plot_list[[1]] <- plot_list[[1]] + theme(
    axis.ticks.y = element_line(),
    axis.line.y = element_line(),
    axis.text.y = element_text(),
    axis.title.y= element_text(),
) + ylab("UMAP_2")
combined_plot <- Reduce(`+`, plot_list)
combined_plot
