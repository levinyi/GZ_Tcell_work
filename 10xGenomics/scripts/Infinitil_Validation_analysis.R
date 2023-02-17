# After screening, print spot count of each screened clonotypes to UMAP.
# input file must contain two columns which are "renamed_clone" and "spots".
# output is a umap plot.

# load the rds.
library(Seurat)
library(ggplot2)
library(tidyverse)

#' create a branded color scale that we can apply to a ggplot2 plot with + scale_color_branded()
# https://www.garrickadenbuie.com/blog/custom-discrete-color-scales-for-ggplot2/
ROOTPATH_COLOURS <- list("black" = "#3E3D57", "violet" = "#FF00F1",
                         "orange" = "#FF996A", "yellow" = "#FFEC00",
                         "green" = "#00F300", "navy" = "#3565AA",
                         "blue" = "#0069FF", "purple" = "#7D2DFD",
                         "red" = "#FC5B68", "skyblue" = "#5FB3D3")
# define some manual color scales
RP_SS_JET <- c(ROOTPATH_COLOURS["black"], ROOTPATH_COLOURS["navy"],
               ROOTPATH_COLOURS["yellow"], ROOTPATH_COLOURS["orange"],
               ROOTPATH_COLOURS["red"], ROOTPATH_COLOURS["violet"])
RP_DS_JET <- c(ROOTPATH_COLOURS["navy"], ROOTPATH_COLOURS["blue"],
               ROOTPATH_COLOURS["skyblue"], ROOTPATH_COLOURS["yellow"],
               ROOTPATH_COLOURS["orange"], ROOTPATH_COLOURS["red"])
RP_GRADIENT <- c(ROOTPATH_COLOURS["navy"],
                 ROOTPATH_COLOURS["yellow"],
                 ROOTPATH_COLOURS["violet"])

#' Create a palette function for branded_colors.
#' what we need from a palette function is a function that takes a single argument n and returns n colors.
rootpath_palette <- function(primary="violet", other="black", reverse=FALSE)
{
  color_names <- names(ROOTPATH_COLOURS)
  # invalid name provided
  stopifnot(primary %in% color_names)

  pal_fun <- function(n_colors) {
    # too many colours
    if (n_colors > length(ROOTPATH_COLOURS)) {
      warning("Requested number of colours exceeds the brand color palette")
    }
    # special 2-color case
    if (n_colors == 2) {
      other <- if (!other %in% color_names) {other}
               else {ROOTPATH_COLOURS[[other]]}
      color_list <- c(other, ROOTPATH_COLOURS[[primary]])
    }
    # generic case
    else {
      start_idx <- which(names(ROOTPATH_COLOURS) == primary)
      if (!reverse) {
        color_list <- c(ROOTPATH_COLOURS[start_idx:length(ROOTPATH_COLOURS)],
                        ROOTPATH_COLOURS[1:start_idx]
                      ) %>% unlist %>% unname
        color_list <- color_list[1:n_colors]
      }
      else {
        color_list <- c(ROOTPATH_COLOURS[start_idx:1],
                        ROOTPATH_COLOURS[length(ROOTPATH_COLOURS):start_idx]
                      ) %>% unlist %>% unname
        color_list <- color_list[1:n_colors]
      }
    }

    return(color_list)
  } # end pal_fun

  return(pal_fun)
} # end rootpath_palette


#' Discrete color gradient with RootPath brand colors, for use in plotting
#' Has same input as rootpath_palette, since it is simply a wrapper that
#' connects it to ggplot2::discrete_scale
scale_colour_rootpath <- function(primary="violet",
                                 other="black",
                                 direction=1,
                                 ...)
{
  ggplot2::discrete_scale(
    aesthetics = "colour",
    scale_name = "RootPath",
    palette = rootpath_palette(primary = primary,
                               other = other,
                               reverse = direction < 0),
    ...
  )
}

#' Convenience function making a multi-color jet gradient over a specified
#' value range. Default colors are a blue-red hue range of 6 brand colours
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


plot_clones <- function(srt_obj, save_path,dim_reduction="umap",
			name="definition"
			){
	plot <- Seurat::DimPlot(object= srt_obj,
				reduction = dim_reduction,
				group.by = name,
				pt.size = 0.2
				) +
			ggtitle("tumor-reactive clones") +
			scale_colour_rootpath()
	if (!is.null(save_path)){
	    ggsave(filename=save_path, plot=plot, device="png")
	}
}



args = commandArgs(T)
srt_obj = readRDS(args[1]) # integrated_seurat.rds

# read spots file
input_file <- read.table(args[2], header=T) # validation file

# Add validation to meta.data
# Add spots to meta.data
clonotypes <- srt_obj[["clonotype_id"]] %>% tibble::rownames_to_column("barcode")
clonotypes$definition <- input_file$definition[match(clonotypes$clonotype_id, input_file$renamed_clone)]
clonotypes <- clonotypes %>% select(c("barcode","definition")) %>% tibble::column_to_rownames("barcode")
srt_obj_temp <- Seurat::AddMetaData(object = srt_obj, metadata = clonotypes, col.name = "definition")
# print(srt_obj_temp@meta.data$definition)
# plot

data = FetchData(srt_obj_temp,vars=c("UMAP_1","UMAP_2","definition"))
data$definition <- paste(data$definition)
data$order <- 3 # NA
data$order[data$definition=="FALSE"] <- 2
data$order[data$definition=="TRUE"] <- 1
print(data[1:5,])
data = data %>% arrange(desc(order))
print(data[1:5,])

plot <- ggplot(data,aes(x=UMAP_1,y=UMAP_2, color=definition))  + 
	geom_point(size=2.5, alpha=0.9)+
	scale_color_manual(values=c("FALSE"="#3D3D55","TRUE"="#FF00F1","NA"="grey"),breaks=c("TRUE","FALSE","NA"))+
	ggtitle("tumor-reactive clones") +
	theme_classic()+
	theme(
	      legend.title=element_blank(),

	)

plot2 <- ggplot(data,aes(UMAP_1,UMAP_2))  + 
 	geom_point(data=filter(data, definition=="NA"),size=2.5, alpha=0.9, color="grey") +
 	geom_point(data=filter(data, definition=="FALSE"),size=2.5, alpha=0.9,color="#3D3D55")+
	geom_point(data=filter(data, definition=="TRUE"),size=2.5, alpha=0.9, color='#FF00F1')+
 	#scale_color_manual(values=c("#3D3D55","#FF00F1","#D3D3D3"))+
 	ggtitle("tumor-reactive clones") +
 	theme_classic() + theme(legend.title=element_blank()
	)

# args[3] is output prefix name eg. SA511-TIL01-10000TCR
save_path <- paste(args[3],"Validation_anlaysis_DimPlot.png",sep="_")
save_path2 <- paste(args[3],"Validation_anlaysis_DimPlot.test.png",sep="_")
# plot_clones(srt_obj_temp, save_path=save_path, name="definition")

ggsave(filename=save_path, plot=plot, device="png")
ggsave(filename=save_path2, plot=plot2, device="png")
