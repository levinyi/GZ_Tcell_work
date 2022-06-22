# A small set of functions to make plotting "standard" plots
# from Seurat objects (in the context of our specific analysis) more
# convenient and less verbose

# Author: Seon Kinrot (RootPath Genomics)

# dependancies
library(Seurat)
library(ggplot2)
library(patchwork)
library(dplyr)

# define company brand palette; based on this example:
# www.garrickadenbuie.com/blog/custom-discrete-color-scales-for-ggplot2/
ROOTPATH_COLOURS <- list("black" = "#3E3D57", "violet" = "#FF00F1",
                         "orange" = "#FF996A", "yellow" = "#FFEC00",
                         "green" = "#00F300", "navy" = "#3565AA",
                         "blue" = "#0069FF", "purple" = "#7D2DFD",
                         "red" = "#FC5B68", "skyblue" = "#5FB3D3")

#' Palette function for RootPath brand
#'
#' Inputs:
#' @param primary - the primary color, used to highlight results
#' @param other - color of secondary "uninteresting" category
#' @param reverse - whether to reverse the order of colours in the palette

#' Output:
#' @param pal_fun - a function that takes a number (n_colors) and returns list
#' of colors of length n_colors, generated from ROOPATH_COLOURS
#'
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
      color_list <- unname(unlist(ROOTPATH_COLOURS[1:n_colors]))
    }

    if (reverse) {
      color_list <- rev(color_list)
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
  discrete_scale(
    aesthetics = "colour",
    scale_name = "RootPath",
    palette = rootpath_palette(primary = primary,
                               other = other,
                               reverse = direction < 0),
    ...
  )
}

#' Function for generating individual feature plots for genes of interest.
#'
#'Inputs:
#' @param gex_seurat - the Seurat object for which plots are made
#' @param plot_dir - the directory in which plots are saved
#' @param name_id - an identifier appearing in all savde file names
#' @param user_genes - a user-supplied list of genes to be plotted. There is an
#' additional default gene list, separated into "DRC" and "Yost" (see
#' first p[ortion of function code). "user_genes" is treated either as
#' additional genes, or as the entire list, depending on the value of
#' "ignore_defaults"
#' @param user_descriptions - a list of "descriptions" with length matching
#' user_genes. Each description associates a gene with a phenotype.
#' @param user_modules - a list of lists. Each sublist is a list of genes
#' whose average expression is treated as a module. This list should be named,
#' with the names being used for metadata and display in the plot.
#' @param um_descriptions - has the same relation to user_modules as
#' user_descriptions is to user_genes.
#' @param axes - the access key for the axes to be plotted. This should be
#' directly callable from the Seurat object's top layer (e.g.
#' gex_seurat[["UMAP_1"]])
#' @para, axes_names - plot labels for the x and y axes, respectively
#' @param colors - a named list specifying the color blend used for expression
#' levels
#' @param ignore_defaults - if true, ignores "DRC" and "Yost" variable and
#' instead plots only user-supplied genes
#' @param print_ - if TRUE, print the plots to the console (in addition to
#' saving to disk)
#' @param w,h,dev,a - the width (inches), height (inches), image format (e.g.
#' "png", "svg") and opacity (0<=alpha<=1) of the image. This is used as input
#' to ggsave with default DPI (300)
#'
#' Returns:
#' @param plot_dir - the directory plots are saved to
#'
gen_feature_plots <- function(gex_seurat, plot_dir,
                              name_id = "FeaturePlots_GexUMAP",
                              user_genes = list(),
                              user_descriptions = list(),
                              user_modules = list(),
                              um_descriptions = list(),
                              axes = c("UMAP_1", "UMAP_2"),
                              axes_names = c("GEX_UMAP_1", "GEX_UMAP_2"),
                              colors = list(
                                          "low" = "lightgrey", #"low" = ROOTPATH_COLOURS[["green"]],
                                          "high" = "blue" #"high" = ROOTPATH_COLOURS[["violet"]]
                                        ),
                              ignore_defaults = FALSE,
                              print_ = FALSE, w=6, h=6, dev="png", a=0.7) {
  # default lists
  drc_genes = c("CCR7", "GZMK", "GZMB", "FASLG", "GZMA", "TCF7", "IL7R",
                "TOX", "KLF2", "LAYN", "ITGAE", "IL2", "TNF", "LAG3",
                "HAVCR2", "SELL", "ICOS", "CD28", "ENTPD1", "IRF2", "IRF3",
                "IRF5", "IRF7", "IRF9", "SOCS1", "STAT1", "STAT2", "TBX21",
                "SCYL3", "CASP3", "CX3CR1", "KLRK1")
  drc_descriptions = c("Naive/ CM", "TE", "Effector/ TE", "TE", "TE", "CM",
                       "Stem/ CM", "Exhaustion?", "Naive/ CM",
                       "Exhaustion in TIL", "TR TIL?", "Effectors",
                       "Effectors", "Effectors", "Effectors", "Naive/ CM",
                       "Co-stim", "Co-stim", "TE/ TR TIL?", "IFN Signaling",
                       "IFN Signaling", "IFN Signaling", "IFN Signaling",
                       "IFN Signaling", "IFN Signaling", "IFN Signaling",
                       "IFN Signaling", "Effectors", "Effectors", "Effectors",
                       "Effectors", "TR TIL?")

  yost_genes = c("CD8A", "CD4", "FOXP3", "EOMES", "KLRD1", "IFNG", "PDCD1")
  yost_descriptions = c("Cytotoxic", "Helper", "T-reg", "Naive/ CM", "TE",
                        "IFN Signaling", "TE")

  modules = list("Yost_exhaustion_score" = c("CTSW", "GNLY", "PRF1", "SLA2",
                                             "HAVCR2", "GZMB", "ITGAE",
                                             "GALNT2", "ACP5", "CXCR6",
                                             "ENTPD1", "KRT86", "TIGIT",
                                             "LAYN", "JAML", "AC092580.4",
                                             "AHI1", "ALOX5AP", "SYNGR2",
                                             "CARD16", "FKBP1A", "GOLIM4",
                                             "SNX9", "TNFRSF18", "ZBED2",
                                             "TNFRSF9", "VCAM1", "CXCL13",
                                             "GEM"),
                  "T_cell_score" = c("CD3D", "CD3E", "CD3G", "CD247"))
  module_descriptions = list("Yost_exhaustion_score" = "Exhaustion",
                             "T_cell_score" = "Tcell")


  # combine feature lists
  if (ignore_defaults) {
    gene_list <- user_genes
    descriptions <- user_descriptions
    modules <- user_modules
    module_descriptions <- um_descriptions
  } else {
    gene_list <- c(drc_genes, yost_genes, user_genes)
    descriptions <- c(drc_descriptions, yost_descriptions,
                      user_descriptions)
    modules <- c(modules, user_modules)
    module_descriptions <- c(module_descriptions, um_descriptions)
  }
  gene_list <- unlist(gene_list)
  descriptions <- unlist(descriptions)

  # add module scores
  temp_seurat <- gex_seurat
  for (module_name in names(modules)) {
    module_feats <- modules[[module_name]]
    retained_feats <- module_feats[module_feats %in% rownames(temp_seurat)]
    if (length(retained_feats) > 1) { # no use in a module of one gene...
      temp_seurat <- AddModuleScore(temp_seurat,
                                    features = retained_feats,
                                    name = module_name)
      gene_list <- c(gene_list, paste0(module_name, "1"))
      descriptions <- c(descriptions,
                        module_descriptions[[module_name]])
    }
  }

  if (!dir.exists(plot_dir)) {
    dir.create(plot_dir)
  }
  # fetch relevant data and plot
  dat <- FetchData(object=temp_seurat, vars=c(axes, gene_list))
  for (ig in seq_len(length.out = length(gene_list))) {
      description <- descriptions[ig]
      gene <- gene_list[ig]
      if (!(gene %in% colnames(dat))) {
          next
      }
      desc_name <- gsub("[/?]", "", gsub(" ", "_", description))
      gplt <- ggplot(dat, mapping = aes(x=.data[[axes[1]]],
                                        y=.data[[axes[2]]])) +
              geom_point(mapping = aes(color=.data[[gene]]), alpha=a) +
              scale_color_gradient(low = colors[["low"]],
                                   high = colors[["high"]]) +
              ggtitle(paste(description, gene)) +
              xlab(axes_names[1]) +
              ylab(axes_names[2])
      if (print_) {
        print(gplt)
      }
      path <- paste0(plot_dir, "/", name_id, "_",
                     desc_name, "_", gene, ".", dev)
      ggsave(filename=path, plot=gplt, width = w, height = h, device = dev)
  }
  return(plot_dir)
}

#' Shows the TCR clonality for all cells in the GEX UMAP.
#' In principle, this can be used for a any purely numerical metadata.
#' Note that this function does not permanently modify the input Seurat
#' object in any way.
#'
#' Inputs:
#' @param gex_seurat - a Seurat object from which to plot
#' clonality_meta - name of the metadata entry in which clonality information
#' is encoded
#' @param title - the plot title
#' @param colors - a named list with entries "low", "mid" and "high", defining
#' the color blend to be used in the plot. The midpoint is identified by the
#' median of uniquely observed clonalities
#' @param alpha - the degree of opacity of each point in the plot
#' @param axes - the access key for the axes to be plotted. This should be
#' directly callable from the Seurat object's top layer (e.g.
#' gex_seurat[["UMAP_1"]])
#' @param axes_names - labels for the x and y axes, respectively
#' @param save_path - if not NULL, specifies a path to save the plot
#' @param w, h, dev - the width and height of the saved image (in inches), as
#' well as the "device" (e.g. "png", "svg"...). This is used as input to
#' ggsave with default DPI (300)
#'
#' Returns:
#' @param plot - a ggplot object of this plot
#'
plot_clonality <- function(gex_seurat, clonality_meta = "TCR_clonality",
                           title = "TCR clonality",
                           colors = list("low"=ROOTPATH_COLOURS[["blue"]],
                                         "mid"=ROOTPATH_COLOURS[["orange"]],
                                         "high"=ROOTPATH_COLOURS[["red"]]),
                           alpha = 0.5,
                           axes = c("UMAP_1", "UMAP_2"),
                           axes_names = c("GEX_UMAP_1", "GEX_UMAP_2"),
                           save_path = NULL, w = 8, h = 6, dev = "png") {
  # fetch data and generate plot
  dat <- FetchData(object=gex_seurat, vars=c(axes, clonality_meta))
  dat[[clonality_meta]][which(is.na(dat[[clonality_meta]]))] <- 0
  midp <- 0.5 * max(dat[[clonality_meta]], na.rm=TRUE)
  plot <- ggplot(dat, mapping = aes(x=.data[[axes[1]]],
                                    y=.data[[axes[2]]])) +
          geom_point(mapping = aes(color=.data[[clonality_meta]],
                                   size=log(.data[[clonality_meta]]+1)),
                     alpha = alpha) +
          scale_color_gradient2(low=colors[["low"]],
                                mid=colors[["mid"]],
                                high=colors[["high"]],
                                midpoint=midp) +
          xlab(axes_names[1]) + ylab(axes_names[2]) + ggtitle(title)
  # save and return
  if (!is.null(save_path)) {
    ggsave(filename=save_path, plot=plot, width=w, height=h, device=dev)
  }
  return(plot)
}

#' Shows the cells from a selected list of clonotypes in the GEX UMAP.
#' In principle, this can be used for a specified selection based on any
#' metadata.
#' Note that this function does not permanently modify the input Seurat
#' object in any way.
#'
#' Inputs:
#' @param gex_seurat - a Seurat object from which to plot
#' @param clones - a clone or list of clone names to be plotted
#' @param clone_metadata - name of the metadata entry in which "clones" can be
#' found
#' @param highlights - a (optional) named list of additional features to be
#' displayed. Each entry needs to be named by the label it will have in the
#' plot and contain a vector with 2 entries: the name of a metadata
#' annotation in gex_seurat and a list or vector of values for that metadata.
#'
#' Example entry:
#' highlights <- list("CD8" = list("seurat_clusters", c(1,3,5,8,9)),
#'                    "CD4" = list("seurat_clusters", c(0,4,6,7,11)),
#'                    "Non-T" = list("seurat_clusters", c(2,10,12)))
#' @param default_val - the value assigned to any cell not selected by clones
#' or highlights
#' @param name - the name of the new metadata entry used for this plot
#' @param point_sz - a vector of length 2: the point size for non-selected and
#' for selected clones
#' @param title - the plot title
#' @param dim_reduction - the name of the dimenstionality reduction to use
#' for the plot (e.g. "umap", "pca", "csp_umap"...)
#' @param axes_names - labels for the x and y axes, respectively
#' @param save_path - if not NULL, specifies a path to save the plot
#' @param w, h, dev - the width and height of the saved image (in inches), as
#' well as the "device" (e.g. "png", "svg"...). This is used as input to
#' ggsave with default DPI (300)
#'
#' Returns:
#' @param plot - a ggplot object of this plot, generated via Seurat DimPlot
#'
plot_clones <- function(gex_seurat, clones, clone_metadata="clonotype_id",
                        highlights=list(), default_val="Other",
                        name="Selected_clones",
                        point_sz=c(0.2, 2.5), title="Selected TCRs",
                        dim_reduction="umap",
                        axes_names=c("GEX_UMAP_1", "GEX_UMAP_2"),
                        save_path=NULL, w=8, h=6, dev="png") {
  # generate metdata
  clonotypes <- gex_seurat[[clone_metadata]][[1]]
  clone_meta <- clonotypes %in% clones
  clone_positions <- which(clone_meta==TRUE)
  clone_meta[clone_positions] <- clonotypes[clone_positions]
  for (feature in names(highlights)) {
    meta_name <- highlights[[feature]][[1]]
    meta_values <- highlights[[feature]][[2]]
    cells <- which(gex_seurat[[meta_name]][[1]] %in% meta_values &
                   clone_meta == FALSE)
    clone_meta[cells] <- feature
  }
  clone_meta[which(clone_meta==FALSE)] <- default_val

  temp_seurat <- AddMetaData(object = gex_seurat,
                             metadata = clone_meta,
                             col.name = name)

  # plot
  ptsz <- rep(point_sz[1], each = length(clone_meta))
  ptsz[clone_positions] <- point_sz[2]
  plot <- DimPlot(object = temp_seurat, reduction = dim_reduction,
                  group.by = name, pt.size = ptsz, order=clones) +
          ggtitle(title) + xlab(axes_names[1]) + ylab(axes_names[2]) +
          scale_colour_rootpath()
  if (!is.null(save_path)) {
    ggsave(filename=save_path, plot=plot,  width=w, height=h, device=dev)
  }
  return(plot)
}

#' Function to plot the distribution of a numerical metadata entry as a violin
#' plot, with some additional visualization customization support
#'
#' Inputs:
#' Returns:
#'
#'
plot_metric_violin <- function(srt_obj, metric, group_by=NULL, save_path=NULL,
                               title="", y_label="",
                               vln_fill=ROOTPATH_COLOURS[["skyblue"]],
                               jitter=0.25, x_lim=c(-1,1), text_size=6,
                               a=0.3, pt_size=1, w=4, h=3, dev="png") {
  dat <- FetchData(srt_obj, vars = c(metric, group_by))
  if (!is.null(group_by)) {
    dat <- dat[order(dat[[group_by]]), ]
  }
  vln_plot <- ggplot(dat, mapping = aes(x = 0, y = .data[[metric]])) +
              geom_violin(fill = vln_fill) +
              geom_jitter(aes(colour = group_by),
                          width = jitter, alpha = a, size = pt_size) +
              scale_colour_rootpath() +
              xlim(x_lim[1], x_lim[2]) +
              theme_bw() +
              ggtitle(title) + ylab(y_label) +
              theme(text = element_text(size = text_size),
                    axis.title.x = element_blank(),
                    axis.text.x = element_blank(),
                    axis.ticks.x = element_blank(),
                    panel.grid.minor = element_blank(),
                    panel.grid.major.x = element_blank())
  if (!is.null(save_path)) {
    ggsave(filename=save_path, plot=vln_plot,
           width = w, height = h, device = dev)
  }
  return(vln_plot)
}
