# A set of wrapper functions that perform specific functions in (post-count)
# single-cell sequencing analysis, using the R package Seurat.
# Most functions take some setting parameters as input, and perform a set of
# operations meant to achieve a "step" in the analysis (e.g. load the data,
# integrate multiple datasets, merge data of different types etc.)

# Author: Seon Kinrot (RootPath Genomics)

# dependancies
library(Seurat)
library(ggplot2)
library(patchwork)
library(dplyr)

#' wrapper to load 10X RNA-seq and CITE-seq data in a single line, while still
#' allowing format flexibility. By default, tries using the folder with the
#' tsv files; if this fails, uses the HDF5 file
#'
#' Inputs:
#' @param prj_dir - the directory containing the filtered_feature_bc_matrix
#'
#' Returns:
#' @param raw_data - the output of one of Seurat's Read10X functions. This will
#' be either a sparse matrix (e.g. for RNA-seq) or a list of such matrices
#' (e.g. for CITE-seq data)
#'
read_10X <- function(prj_dir) {
  raw_data <- tryCatch(
    {
      # we want the filtered_feature_bc_matrix folder, but sometimes its name
      # contains other things (e.g. project name)
      data_path <- grep(pattern = "filtered_feature_bc_matrix",
                        x = list.dirs(prj_dir, full.names = T, recursive = T),
                        value = T,
                        ignore.case = T)
      Read10X(data.dir = data_path)
    },
    error = function(cond) {
      # try to find the HDF5 file and load it
      hdf_fl <- list.files(path = prj_dir,
                           pattern = "filtered_feature_bc_matrix.h5",
                           full.names = T, recursive = T)
      Read10X_h5(filename = hdf_fl)
    }
  )
  return(raw_data)
}

#' wrapper function that loads single-cell RNA-seq 10x data into a Seurat
#' object, given the folder it resides in
#'
#' Inputs:
#' @param prj_dir - the directory where your data is. Should be the parent
#' directory of the filtered_feature_bc_matrix folder. Also please note that
#' Seurat needs to have the *compressed* files in that sub-folder.
#' @param prj_nm - optional name for your project to save in Seurat metadata.
#' If NULL (default) it will give the project the basename of the folder
#' @param min_freq, min_min_cells, max_min_cells - these variables are used to
#' determine the minimal number of cells a gene is expressed in to be included
#' in the final gene list. load_seurat will set the minimal number of cells to
#' min.min.cells (default 3), unless the target minimal expression rate
#' (min.freq, default 1e-3) corresponds to 4 or more cells. Then it will
#' enforce that target minimal expression rate, up to a certain maximal number
#' (max.min.cells, default 10). This means even if the target minimum
#' expression rate is higher, we will keep any gene expressed in at least
#' max.min.cells cells.
#' @param min_feats - minimal number of features in order to keep a cell.
#' Default is 0, since this is later filtered more carefully by downstream QC
#' @param rmv_grep, rmv_exact - lists of genes to be removed from further
#' analysis, irrespective of expression levels. rmv.grep are expressions used
#' to identify groups of genes containing them (e.g., "TRAV" to remove all
#' TCR-alpha chain variable genes). rmv.exact is a list of genes requiring an
#' exact match in gene names (e.g., a list of confounding genes. Note that this
#' removes genes from the Seurat object's count list, and does not regress
#' other effects associated with their expression)
#' @param fix_bcodes - whether or not to remove the "-1" at the end of each
#' barcode
#' @param verbose - if True, the function will print reports on number of genes
#' removed at each step, etc
#
#' Output:
#' @param srt_obj - a Seurat object ready for QC
#'
load_seurat <- function(prj_dir, prj_nm=NULL, min_freq=1e-3, min_min_cells=3,
                        max_min_cells=10, min_feats=0, rmv_grep=NULL,
                        rmv_exact=NULL, fix_bcodes=TRUE, verbose=FALSE) {
  # make sure this is loading a single dataset
  if (length(prj_dir)!=1) {
    cat("Did not provide a single folder!\n")
    return(NULL)
  }
  if (!dir.exists(prj_dir)) {
    cat(prj_dir, "does not exist!\n")
    return(NULL)
  }

  raw_data <- read_10X(prj_dir)
  # if dim is NULL, assumes this is CITE-seq data, and is a list of length 2
  adt_raw <- NULL
  if (is.null(dim(raw_data))) {
    adt_raw <- raw_data[["Antibody Capture"]]
    raw_data <- raw_data[["Gene Expression"]]
  }

  if (verbose) {
    cat("Loaded raw data, with", raw_data@Dim[2], "cells and", raw_data@Dim[1],
        "genes\n", sep = " ")
  }

  # remove unwanted genes; this behavior is deprectated, and feature removal is
  # now implemented in the variable gene step. Retained for compatibility
  keep_genes <- filter_genes(genes = rownames(raw_data),
                             rmv_grep = rmv_grep,
                             rmv_exact = rmv_exact)
  filtered_data <- raw_data[keep_genes, ]

  # fix barcodes
  if (fix_bcodes) {
    colnames(filtered_data) <- gsub(pattern = "-1", replacement = "",
                                    x = colnames(filtered_data))
    if (!is.null(adt_raw)) {
      colnames(adt_raw) <- gsub(pattern = "-1", replacement = "",
                                x = colnames(adt_raw))
    }
  }
  if (verbose) {
    cat("Unwanted genes filtered; removed", sum(!keep_genes),
        "genes overall\n", sep = " ")
  }

  #determine minimal number of expressing cells for a gene to be retained
  ncells <- filtered_data@Dim[2]
  min_cells <- min(
    max(floor(ncells * min_freq), min_min_cells),
    max_min_cells
  )
  # generate Seurat object
  if (is.null(prj_nm)) {
    prj_nm <- basename(prj_dir)
  }
  srt_obj <- CreateSeuratObject(counts = filtered_data, project = prj_nm,
                                min.cells = min_cells, min.features = min_feats)
  if (!is.null(adt_raw)) {
    adt_assay <- CreateAssayObject(counts = adt_raw)
    srt_obj[["ADT"]] <- adt_assay
  }

  if (verbose) {
    cat("Generated Seurat object (min.cells =", min_cells,
        "and min.features = ", min_feats,
        "). Final matrix size (Genes x cells):",
        srt_obj@assays$RNA@counts@Dim[1], "x", srt_obj@assays$RNA@data@Dim[2],
        "\n", sep = " ")
  }

  return(srt_obj)
}

#' A helper function for removing genes from a list of gene names.
#'
#' Inputs:
#' @param genes - a list of gene names
#' @param rmv_grep - a list of regular expressions. Any gene matching these
#' will be discarded
#' @param rmv_exact - a list of gene names to be removed if matched exactly
#' Returns:
#' @param keep_genes - a binary vector of the same length as genes. Entries are
#' TRUE/ FALSE depending on whether the gene is to be retained
#'
filter_genes <- function(genes, rmv_grep = NULL, rmv_exact = NULL) {
  keep_genes <- rep(TRUE, each = length(genes))
  for (gene in rmv_grep) {
    keep_genes[grep(gene, genes, value = F)] <- FALSE
  }
  for (gene in rmv_exact) {
    keep_genes[which(genes == gene)] <- FALSE
  }
  return(keep_genes)
}


#' runs standard QC based on number of UMIs, number of unique genes and
#' percentage of mitochondrial genes present in each cell.
#' This function is meant to be semi-interactive/ iterative. You can run it
#' first just to plot the distributions of QC parameters, then repeatedly run
#' with harsher QC cutoffs until you are satisfied. Since it only returns the
#' list of cells that passed, enforcing a subset is up to the user, and done
#' externally to this function.
#'
#' Inputs:
#' @param srt_obj - a Seurat object to undergo QC
#' @param feature_sets - a named list of feature sets to be considered for QC.
#' The namee are used for metadata entry, and the content of each list entry is
#' a pattern to be used for Seurat's PercenTageFeatureSet function (e.g. ^MT-).
#' Default is to include Mitochondrial and Ribosomal gene content.
#' @param min_fsets, max_fsets - each of these is a vector of minimum (or
#' maximum) values allowed for each of the feature sets in feature_sets.
#' @param min_feats, max_feats etc. - cutoff values for QC. Names should be
#' self-explanatory.
#' @param species - the species being sequenced. Changes the pattern search for
#' Mitochondrial and Ribosomal genes, as well as cell cyel scoring. Not meant
#' to be used in conjkuction with custom feature sets.
#' @param plot_ - a logical value input. If TRUE, the function will genrate the
#' following plots: 1) scatter plot of UMI vs. Genes (with Mitochondrial
#' content as color scale) 2) scatter plot of UMI vs. Mitochondrial content
#' 3) histogram for gene counts per cell and 4) histogram of Mitochondrial
#' content per cell. All plots are for the population of post-subset cells.
#' @param return_plot - a logical value. If TRUE (and plot_ is also TRUE), the
#' function will return both a Seurat object and the QC plot
#' @param return_object - a logical value. If TRUE returns the Seurat object,
#' in addition to any other returned values
#' @param verbose - wehether to print QC results (how many cells are retained)
#' @param score_cc - whether to add cell cycle scoring to metadata. This uses
#' Seurat's CellCycleScoring function, and does not specify whether or in what
#' form to regress cell cycle effects. Only relevant if return_object is TRUE
#'
#' Output:
#' @param ret_list - a list of returned objects containing a subset of the
#' following
#' @param keep_cells - a list of cells that passed the QC thresholds specified.
#' Access under name "keep_barcodes"
#' @param combined_plot - if both plot_ and return.plot are TRUE, it will also
#' return the QC plot as a ggplot object. Access under "plot"
#' @param srt_new - the input srt_obj, with added Mitochondrial gene percentage
#' and filtered by the input parameters. Returned only if return_object is TRUE.
#' Access under "filtered_seurat"
#'
run_qc <- function(srt_obj,
                   feature_sets=list("pct_mito" = "^MT-",
                                     "pct_ribo" = "^RP[SL]"),
                   min_fsets=c(0,0), max_fsets=c(100, 100),
                   min_feats=0, max_feats=Inf,
                   min_umi=0, max_umi=Inf,
                   species='human',
                   plot_=TRUE, return_plot=FALSE, return_object=FALSE,
                   score_cc=TRUE, verbose=FALSE) {
  # start fetching the data for QC
  srt_new <- srt_obj
  feats_ <- FetchData(object = srt_new, vars = "nFeature_RNA")
  umi_ <- FetchData(object = srt_new, vars = "nCount_RNA")
  # specifying mouse as species overrides manual feature_sets
  if (species == 'mouse') {
    feature_sets = list("pct_mito" = "^mt-",
                        "pct_ribo" = "^Mrp[sl]")
  }
  qc_sets <- list()
  # add feature sets for QC
  for (feature_name in names(feature_sets)) {
    srt_new <- PercentageFeatureSet(object = srt_new,
                                    pattern = feature_sets[[feature_name]],
                                    col.name = feature_name)
    qc_sets[[feature_name]] <- FetchData(object = srt_new, vars = feature_name)
  }
  # subset the data
  # See https://github.com/satijalab/seurat/issues/1053 for discussion of how
  # to subset data in a programmatic way
  keep_cells <- feats_ <= max_feats & feats_ >= min_feats &
                umi_ <= max_umi & umi_ >= min_umi
  for (ifeature in seq_len(length(qc_sets))) {
    keep_cells <- keep_cells &
                  qc_sets[[ifeature]] <= max_fsets[ifeature] &
                  qc_sets[[ifeature]] >= min_fsets[ifeature]
  }
  keep_cells <- which(keep_cells)
  ret_list <- list("keep_barcodes"=keep_cells)
  srt_new <- srt_new[, keep_cells]

  # print update on retained cells
  if (verbose) {
    cat("Retained", length(keep_cells),
        "cells out of", dim(srt_obj)[2], "\n",
        sep=" ")
  }

  if (return_object) {
    if (score_cc) {
      srt_new <- score_cell_cycle(gex_seurat=srt_new, species=species)
    }
    ret_list <- c(ret_list, list("filtered_seurat"=srt_new))
  }

  if (plot_) {
    plot_res <- plot_qc(srt_new, names(feature_sets))

    if (return_plot) {
      ret_list <- c(ret_list, list("plots"=plot_res[["list"]]))
    }
    else {
      print(plot_res[["plot"]])
    }
  }

  return(ret_list)
}

#' A helper function that generates the QC plot for run_qc.
#' In principle, it would be best to move this function to the PlotFunctions
#' module (and also modify it so it uses the brand colours), but this is
#' somewhat cumbersomean due to path handling issues when this file is sourced
#' from outside the project (e.g. a notebook).
#' Having this function "live" here is an interim solution.
#'
#' Inputs:
#' @param srt_new - a Seurat object containing the data to be plotted
#' @param - feature_sets - the names of metadata to be plotted (in addition to
#' UMI and gene counts)
#' @param a - alpha for the plots
#' @param bins - number of bins in histograms
#' @param mito_grad - the UMI vs. gene abundance plots encodes dot (barcode)
#' color by the first entry in feature_sets (the default is Mitochondrial gene
#' content). This defines the colors to be used for that gradient. These colors
#' are also used for the fill colou of the estimated density of histograms
#' @param hist_colors - color and fill for the histogram bars
#' Returns:
#' @param plot_list - a list of all individual plots made by this function
#' @param combined_plot - a patchwork object containing all QC plots in a grid
#'
plot_qc <- function(srt_new,
                    feature_sets=c("pct_mito", "pct_ribo"),
                    a=0.1, bins=50,
                    mito_grad=c("navy", "red"),
                    hist_colors=c("black", "white")) {
  # fetch QC data
  qc_data <- FetchData(object = srt_new,
                       vars = c("nCount_RNA", "nFeature_RNA"))
  qc_data <- log10(qc_data + 1)
  names(qc_data) <- c("log10_UMI", "log10_genes")
  for (feature_name in feature_sets) {
    column <- FetchData(object = srt_new, vars = feature_name)[, 1]
    qc_data[[feature_name]] <- column
  }
  plot_list <- list()
  # plot umi vs. gene plot and start to compile composite QC plot
  color_feat <- feature_sets[1]
  umi_gene_plot <- ggplot(data = qc_data,
                          mapping = aes(x = log10_genes, y = log10_UMI)) +
                   geom_point(alpha = a,
                              mapping = aes(colour = .data[[color_feat]])) +
                   scale_color_gradient(low = mito_grad[1],
                                        high = mito_grad[2])
  plot_list[["umi_gene_scatter"]] <- umi_gene_plot
  combined_plot <- umi_gene_plot
  # add umi vs. feature (e.g. Mitochondrial genes) plots
  for (feature_name in feature_sets) {
    umi_feat_plot <- ggplot(data = qc_data,
                            mapping = aes(x = .data[[feature_name]],
                                          y = log10_UMI)) +
                     geom_point(alpha = a)
    plot_name <- paste("umi", feature_name, "scatter", sep="_")
    plot_list[[plot_name]] <- umi_feat_plot
    combined_plot <- combined_plot + umi_feat_plot
  }
  # start adding the histograms
  gene_hist <- ggplot(data = qc_data, mapping = aes(x = log10_genes)) +
               geom_histogram(bins = bins, mapping = aes(y = ..density..),
                              colour = hist_colors[1], fill = hist_colors[2]) +
               geom_density(alpha = a, fill = mito_grad[2])
  plot_list[["gene_hist"]] <- gene_hist
  combined_plot <- combined_plot + gene_hist
  # add histograms for each feature set
  for (feature_name in feature_sets) {
    feat_hist <- ggplot(data = qc_data,
                        mapping = aes(x = .data[[feature_name]])) +
                 geom_histogram(bins = bins,
                                mapping = aes(y = ..density..),
                                colour = hist_colors[1],
                                fill = hist_colors[2]) +
                 geom_density(alpha = a, fill = mito_grad[1])
    plot_name <- paste(feature_name, "hist", sep="_")
    plot_list[[plot_name]] <- feat_hist
    combined_plot <- combined_plot + feat_hist
  }
  return(list("list"=plot_list, "plot"=combined_plot))
}

#' A master function to perform integration on a list of GEX 10X datasets. This
#' function supports multiple methods of normalization and integration, and
#' performs all steps from renaming barcodes to manipulations on individual
#' datesets, anchor finding and complete data integration.
#' Inputs:
#' @param srt_list - a list of Seurat objects, one for each experiment to be
#' integrated. Can be named or unnamed.
#' @param normalization_method - specifies the normalization to perform on
#' individual datasets (e.g. "LogNormalize", "SCT", "CLR")
#' @param integration_method - the method used to calculate anchors for
#' integration (e.g. "cca", "rpca")
#' @param k_anchors - the number of anchors to be calculated for local
#' correction (see Seurat documentation). May be modified in some cases,
#' especially when using rpca
#' @param regress - a list of metadata entries to regress during scaling
#' @param cc_regress - whether and in what form to regress cell cycle effects.
#' Can take on the values NULL, "diff", or "all" (see
#' normalize_and_find_var_feats documentation)
#' @param nfeatures - the number of features to be used for integration. If
#' NULL (default), Seurat's default for the chosen normalization is used
#' @param rmv_grep, rmv_exact - genes to exclude from variable feature list
#' @param ref_data - specify which datasets to use as a reference for
#' correction. If NULL (default), no reference is used (and all pairwise
#' comparisons) are made
#' @param readd_counts - specifyies whether to add information to the "count"
#' slot of the "integrated" assay generated
#'
#' Output:
#' @param integrated_seurat - a single Seurat object containing the final,
#' integrated data
#'
integrate_datasets <- function(srt_list, normalization_method = "SCT",
                               integration_method = "rpca", k_anchors = 5,
                               regress = NULL, cc_regress = NULL,
                               nfeatures = 3000,
                               rmv_grep = NULL, rmv_exact = NULL,
                               ref_data = NULL, readd_counts = FALSE) {
  # set global future parameters
  options(future.globals.maxSize = 16*2^30)
  options(future.seed = TRUE)
  # rename barcodes
  samples <- names(srt_list)
  if (is.null(samples)) {
   samples <- seq_len(length(samples))
  }
  for (idata in seq_len(length(srt_list))) {
    srt_list[[idata]] <- RenameCells(object = srt_list[[idata]],
                                     add.cell.id = samples[idata])
  }

  # if we re-add the counts, we will need to know what the raw assay was
  # we will assume this is the same for all samples (and not check this)
  raw_assay <- DefaultAssay(srt_list[[1]])
  if (normalization_method == "SCT") {
    raw_assay <- "SCT"
  }

  # normalization and feature selection. Note that the underlying operations
  # performed for SCT and log normalization vary here
  srt_list <- lapply(X = srt_list, FUN = normalize_and_find_var_feats,
                     normalization_method = normalization_method,
                     regress = regress, cc_regress = cc_regress,
                     nfeatures = nfeatures, rmv_grep = rmv_grep,
                     rmv_exact = rmv_exact, do_scale = FALSE)
  # if regressing cell cycle effects and not regression was not applied in the
  # above line (depends on algorithm choice), regress needs to be upadted
  regress <- update_regress(regress=regress, cc_regress=cc_regress)

  # select integration features
  integration_features <- SelectIntegrationFeatures(object.list = srt_list,
                                                    nfeatures = nfeatures,
                                                    verbose = FALSE)
  # remove unwanted genes
  keep_feats <- filter_genes(genes = integration_features,
                             rmv_grep = rmv_grep,
                             rmv_exact = rmv_exact)
  integration_features <- integration_features[keep_feats]

  # additional step for SCT
  if (normalization_method == "SCT") {
    srt_list <- PrepSCTIntegration(object.list = srt_list,
                                   anchor.features = integration_features,
                                   verbose = FALSE)
  }
  # run PCA if necessary
  if (integration_method == "rpca") {
    # need to scale data for log normalizations
    if (normalization_method != "SCT") {
      srt_list <- lapply(X = srt_list, FUN = ScaleData,
                         features = integration_features,
                         vars.to.regress = regress, verbose = FALSE)
    }
    srt_list <- lapply(X = srt_list, FUN = RunPCA,
                       features = integration_features, verbose = FALSE)
  }

  # make sure reference format is correct, if applicable
  if (!is.null(ref_data)) {
    if (type(ref_data[1] == "character")) {
      temp_ <- vector(mode = "list", length = length(ref_data))
      for (idata in seq_len(length(ref_data))) {
        data_name <- ref_data[idata]
        temp_[idata] <- which(ref_data == data_name)
      }
      ref_data <- temp_
      rm(temp_)
    }
  }
  # identify acceptable normalization name ("CLR" is not accepted)
  norm_str <- normalization_method
  if (normalization_method == "CLR") {
    norm_str <- "LogNormalize"
  }
  # find integration anchors and integrate
  integration_anchors <- FindIntegrationAnchors(object.list = srt_list,
                                reference = ref_data,
                                normalization.method = norm_str,
                                anchor.features = integration_features,
                                reduction = integration_method,
                                k.anchor = k_anchors, verbose = FALSE)
  integrated_seurat <- IntegrateData(anchorset = integration_anchors,
                                   normalization.method = norm_str,
                                   verbose = FALSE)

  # perform scaling and regression, if this was not done yet
  if (normalization_method != "SCT" & integration_method != "rpca") {
    integrated_seurat <- ScaleData(object = integrated_seurat,
                                   vars.to.regress = regress, verbose = FALSE)
  }

  if (readd_counts) {
    raw_data <- integrated_seurat[[raw_assay]]@counts
    raw_data <- raw_data[integration_features, ]
    integrated_seurat[["integrated"]]@counts <- raw_data
  }

  return(integrated_seurat)
}


#' converts mouse gene list to their human conterparts
convert_to_human <- function(mouse_genes){
  # add this line here to allow import without it
  library(biomaRt)
  human_db <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
  mouse_db <- useMart("ensembl", dataset = "mmusculus_gene_ensembl")
  mgenes_v2 <- getLDS(attributes = c("mgi_symbol"),
                      filters = "mgi_symbol",
                      values = mouse_genes,
                      mart = mouse_db,
                      attributesL = c("hgnc_symbol"),
                      martL = human_db,
                      uniqueRows = TRUE)
  human_genes <- unique(mgenes_v2[, 2])
  return(human_genes)
}


#' converts human gene list to their mouse conterparts
convert_to_mouse <- function(human_genes){
  # add this line here to allow import without it
  library(biomaRt)
  human_db <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
  mouse_db <- useMart("ensembl", dataset = "mmusculus_gene_ensembl")
  hgenes_v2 <- getLDS(attributes = c("hgnc_symbol"),
                      filters = "hgnc_symbol",
                      values = human_genes,
                      mart = human_db,
                      attributesL = c("mgi_symbol"),
                      martL = mouse_db,
                      uniqueRows = TRUE)
  mouse_genes <- unique(hgenes_v2[, 2])
  return(mouse_genes)
}


#' Simple function that adds cell cycle scores in one line
score_cell_cycle <- function(gex_seurat, species='human') {
  s_genes <- cc.genes.updated.2019$s.genes
  g2m_genes <- cc.genes.updated.2019$g2m.genes
  if (species == 'mouse') {
    s_genes <- convert_to_mouse(s_genes)
    g2m_genes <- convert_to_mouse(g2m_genes)
  }
  srt_new <- CellCycleScoring(object = gex_seurat,
                              s.features = s_genes,
                              g2m.features = g2m_genes)
  return(srt_new)
}


#' helper function for updating the list of regressed genes for cell cycle
#' effects. This assumes scoring has already been done, and simply updates the
#' list of features to regress out.
update_regress <- function(regress, cc_regress) {
  # edge case 1 - no additional regression
  if (is.null(cc_regress)) {
    return(regress)
  }
  # invalid value case
  if (!(cc_regress %in% c("diff", "all"))) {
    cat("Invalid value for cc_regress! Please specify NULL,",
        "'all' or 'diff'\n", sep=" ")
    return(regress)
  }
  # use difference of S-G2M
  if (cc_regress == "diff") {
    cc_feats <- "CC_diff"
  }
  # use S and G2M scores separately
  if (cc_regress == "all") {
    cc_feats <- c("S.Score", "G2M.Score")
  }
  # update list of variables to regress
  if (is.null(regress)) {
    regress <- cc_feats
  }
  else {
    regress <- c(regress, cc_feats)
  }
  return(regress)
}


#' Helper function that performs normalization and finds variable features
#' for a Seurat object, based on some user-specified input
#'
#' Inputs:
#' @param srt_obj - Seurat object to normalize
#' @param normalization_method - method for normalization (e.g. "SCT")
#' @param regress - list of variables to regress
#' @param cc_regress - whether and in what form to regress cell cycle effects.
#' can be one of the following values - NULL means no regression; "all" means
#' to regress all cell cycle effects (both G2M and S scores) and "diff" means
#' regressing only the difference between S and G2M scores
#' @param nfeatures - how many features to use for variable feature selection.
#' If NULL (default), 3000 features are selected for SCT and 2000 for other
#' methods
#' @param rmv_grep, rmv_exact - genes to be specifically excluded from variable
#' gene list (either by regular expression search or exact match)
#' @param do_scale - whether to scale the data, if applyiong log normalization.
#' Necessary if you want to apply regression without SCT
#'
#' Output:
#' @param srt_obj - a normalized and variable-feature-selected Seurat object.
#' If using SCTranform, scaling and regression of confounding factors is also
#' performed
normalize_and_find_var_feats <- function(srt_obj, normalization_method,
                                         regress = NULL, cc_regress = NULL,
                                         nfeatures = NULL,
                                         rmv_grep = NULL, rmv_exact = NULL,
                                         do_scale = TRUE) {
  # incorporate cell-cycle regression choice
  if (!is.null(cc_regress)) {
    # check for valid cc_regress value
    if (!(cc_regress %in% c("all", "diff"))) {
      cat("Invalid value for cc_regress! Please specify NULL,",
          "'all' or 'diff'\n", sep=" ")
      return(NULL)
    }
    # first, make sure cell cycle is scored
    srt_obj <- score_cell_cycle(gex_seurat=srt_obj)
    if (cc_regress == "diff") {
      srt_obj[["CC_diff"]] <- srt_obj$S.Score - srt_obj$G2M.Score
    }
    # update regress
    regress <- update_regress(regress=regress, cc_regress=cc_regress)
  }

  # SCT normalization
  if (normalization_method == "SCT") {
    if (is.null(nfeatures)) {
      nfeatures <- 3000
    }
    srt_obj <- SCTransform(object = srt_obj, vars.to.regress = regress,
                           variable.features.n = nfeatures, verbose = FALSE)
    }
  # Log normalization or CLR
  else {
    if (is.null(nfeatures)) {
      nfeatures <- 2000
    }
    # set margin on which to normalize data (1 is genes, 2 is cells)
    margin <- 1
    if (normalization_method == "CLR") {
      margin <- 2
    }
    # perform normalization
    srt_obj <- NormalizeData(object = srt_obj,
                             normalization.method = normalization_method,
                             margin = margin, verbose = FALSE)
    srt_obj <- FindVariableFeatures(object = srt_obj, selection.method = "vst",
                                    nfeatures = nfeatures, verbose = FALSE)
    if (do_scale) {
      srt_obj <- ScaleData(object = srt_obj, vars.to.regress = regress,
                           verbose = FALSE)
    }
  }

  # filter out genes from variable features
  aa <- srt_obj@active.assay
  var_feats <- srt_obj[[aa]]@var.features
  keep_feats <- filter_genes(genes = var_feats,
                             rmv_grep = rmv_grep,
                             rmv_exact = rmv_exact)
  var_feats <- var_feats[keep_feats]
  srt_obj[[aa]]@var.features <- var_feats

  return(srt_obj)
}


#' Function that takes a post-QC Seurat object and performs clustering, as well
#' as the necessary preparations (e.g. dimensionality reduction, neighbor graph
#' generation etc). This function is a "shortcut", in the sense that it skips
#' usual steps to determine the appropriate number of dimensions with which to
#' perform clustering etc.
#'
#' Inputs:
#' @param srt_obj - a Seurat object on which to run clustering
#' @param pca_dims - the number of PCA dimensions to calculate for the PCA
#' reduction
#' @param cluster_dims - the PCA dimensions on which to perform UMAP
#' calculation, neighbor finding and clustering
#' @param cluster_res - the resolution/s at which to perform graph-based
#' clustering
#'
#' Outputs:
#' @param srt_obj - the original Seurat object with additional dimensionality
#' reductions etc calculated
#'
cluster_seurat <- function(srt_obj, pca_dims = 50, cluster_dims = 1:30,
                           cluster_res = c(0.4, 0.8, 1.2, 1.6, 2)) {
  srt_obj <- RunPCA(object = srt_obj, npcs = pca_dims, verbose = FALSE)
  srt_obj <- RunUMAP(object = srt_obj, dims = cluster_dims, verbose = FALSE)
  srt_obj <- FindNeighbors(object = srt_obj, dims = cluster_dims,
                           verbose = FALSE)
  srt_obj <- FindClusters(object = srt_obj, resolution = cluster_res,
                          verbose = FALSE)
  return(srt_obj)
}


#' Function to add a Seurat object as a new assay to an existing one. Original
#' intended use is for adding CSP data to a GEX Seurat. Naming convention
#' reflects this use case.
#'
#' Inputs:
#' @param main_srt - the Seurat object to which additional data will be added
#' @param secondary_srt - the Seurat object from which data is transferred
#' @param assay - an assay name from secondary_srt to be added to
#' main_srt. These should be names of assays, specifically (e.g.
#' "integrated"), which can be accessed directly from the top
#' layer of main_srt (e.g. main_srt[["integrated"]])
#' @param new_name, new_key - the new name and key for the assay after being
#' added to main_srt. If NULL (default), the project name will be concatenated
#' at the head of the original assay name/ key
#' @param add_reductions - a boolean variable. If TRUE, goes through all
#' dimensionality reduction objects in secondary_srt and adds those that have
#' their @assay.used parameter set to "assay"
#' @param filter_main - whether to keep only cells with data in both objects.
#' If TRUE, only barcodes shared across the two objects are retained. Otherwise,
#' all barcodes in main_srt are kept and the ones missing from secondary_srt
#' assgined NA values in the new assay
#'
#' Output:
#' @obj main_srt - a Seurat object containing all information in the input
#' main_srt, plus the elements copied from secondary_srt
#'
add_assay <- function(main_srt, secondary_srt, assay,
                      new_name=NULL, new_key=NULL,
                      add_reductions=TRUE, filter_main=FALSE) {
  # resolve assay name and key
  if (is.null(new_name)) {
    prj_nm <- secondary_srt@project.name
    new_name <- paste(prj_nm, assay, sep="_")
  }
  if (is.null(new_key)) {
    # keys don't allow underscores, except as last character
    new_key <- paste0(tolower(gsub("_", "", new_name)), "_")
  }
  # assays can only be added when all cells match - identify these cells
  main_bcodes <- colnames(main_srt)
  sec_bcodes <- colnames(secondary_srt)
  common_bcodes <- main_bcodes[which(main_bcodes %in% sec_bcodes)]
  sec_subset <- secondary_srt[, common_bcodes]
  sec_assay <- sec_subset[[assay]]
  # simple case - main_srt is subset, as well
  if (filter_main) {
    sec_assay@key <- new_key
    main_subset <- main_srt[, common_bcodes]
    main_subset[[new_name]] <- sec_assay
  }
  # involved case - we need to add NA for each cell in main_srt not present in
  # secondary_srt
  else {
    main_subset <- main_srt
    dropped_bcodes <- main_bcodes[which(!(main_bcodes %in% sec_bcodes))]
    sec_feats <- rownames(sec_subset)
    sec_bcodes <- c(common_bcodes, dropped_bcodes)
    nrows <- length(sec_feats)
    ncols <- length(sec_bcodes)
    # start to generate new assay
    new_matrix <- matrix(data = rep(NA, each=nrows*ncols), nrow=nrows)
    rownames(new_matrix) <- sec_feats
    colnames(new_matrix) <- sec_bcodes
    # @counts entry
    has_counts <- length(sec_assay@counts) > 0
    if (has_counts) {
      new_counts <- new_matrix
      new_counts[, common_bcodes] <- as.matrix(sec_assay@counts[, common_bcodes])
      new_assay <- CreateAssayObject(counts=new_counts)
    }
    # @data entry
    has_data <- length(sec_assay@data) > 0
    # if no counts or data, something is wrong
    if (!has_counts & !has_data) {
      cat("No counts or data found in specified assay!")
      return(NULL)
    }
    # if no data, nothing needs to be done
    if (has_data) {
      new_data <- new_matrix
      new_data[, common_bcodes] <- as.matrix(sec_assay@data[, common_bcodes])
      # create the assay, or add to the one with the counts
      if (!has_counts) {
        new_assay <- CreateAssayObject(data=new_data)
      }
      else {
        new_assay@data <- new_data
      }
    }
    # @scale.data entry
    has_scale <- length(sec_assay@scale.data) > 0
    if (has_scale) {
      # scale.data may be for a subset of features (e.g. variable features)
      scale_feats <- rownames(sec_assay@scale.data)
      new_scale <- new_matrix[scale_feats, ]
      new_scale[, common_bcodes] <- as.matrix(
                                        sec_assay@scale.data[, common_bcodes])
      new_assay@scale.data <- new_scale
    }
    # non-data entries
    new_assay@key <- new_key
    new_assay@var.features <- sec_assay@var.features
    new_assay@meta.features <- sec_assay@meta.features
    new_assay@misc <- sec_assay@misc
    new_assay@assay.orig <- sec_assay@assay.orig
    # add assay to main_srt
    main_subset[[new_name]] <- new_assay
  }
  # end assay creation scenario

  # dimensionality reduction loop; this is done irrespective of assay creation
  if (add_reductions) {
    sec_reductions <- names(secondary_srt@reductions)
    assays_used <- sapply(X=sec_reductions,
                          FUN=function(red) {
                                secondary_srt[[red]]@assay.used
                              })
    for (reduction in sec_reductions) {
      if (assays_used[reduction] == assay) {
        main_subset <- add_reduction(main_srt=main_subset,
                                     secondary_srt=secondary_srt,
                                     reduction=reduction,
                                     new_name=paste(tolower(new_name),
                                                    reduction,
                                                    sep="_"),
                                     assay_name=new_name,
                                     filter_main=filter_main)
      }
    }
  }

  return(main_subset)
}


#' Function to add a dimensionality reduction from one Seurat object to another.
#' For example, this can be used to add the UMAP or tSNE reduction from the CSP
#' to the GEX portion of a CITE-seq experiment.
#'
#' Inputs:
#' @param main_srt - the Seurat object to which additional data will be added
#' @param secondary_srt - the Seurat object from which data is transferred
#' @param reduct - the object name by which the dimenstional reduction can be
#' accessed in secondary_srt (e.g. "umap")
#' @param new_name, new_key - the new name and key for the assay after being
#' added to main_srt. If NULL (default), the project name will be concatenated
#' at the head of the original assay name/ key
#' @param assay_name - the name of the assay used to compute the dimensional
# reduction. For example, if this assay is being added to main_srt via
#' add_assay, this should be the same as new_name in that function. If NULL,
#' the assay.used is et to the same as the original reduction, with the project
#' name prepended.
#' @param filter_main - whether to keep only cells with data in both objects.
#' If TRUE, only barcodes shared across the two objects are retained. Otherwise,
#' all barcodes in main_srt are kept and the ones missing from secondary_srt
#' assgined NA values in the new assay
#'
#' Output:
#' @obj main_srt - a Seurat object containing all information in the input
#' main_srt, plus the elements copied from secondary_srt
#'
add_reduction <- function(main_srt, secondary_srt, reduction,
                          new_name=NULL, new_key=NULL, assay_name=NULL,
                          filter_main=FALSE) {
  # resolve name and key
  if (is.null(new_name)) {
    prj_nm <- secondary_srt@project.name
    new_name <- paste(prj_nm, reduction, sep="_")
  }
  if (is.null(new_key)) {
    # keys don't allow underscores, except as last character
    new_key <- paste0(toupper(gsub("_", "", new_name)), "_")
  }
  # dimensionality reductions can only be added when all cells match
  main_bcodes <- colnames(main_srt)
  sec_bcodes <- colnames(secondary_srt)
  common_bcodes <- main_bcodes[which(main_bcodes %in% sec_bcodes)]
  sec_subset <- secondary_srt[, common_bcodes]
  sec_reduct <- sec_subset[[reduction]]
  # assay name for this reduction
  if (is.null(assay_name)) {
    assay_name <- paste(sec_subset@project.name, sec_reduct@assay.used, sep="_")
  }
  # simple case - main_srt is subset, as well
  if (filter_main) {
    sec_reduct@key <- new_key
    # change dimension names to avoid ambiguity
    ncols <- dim(sec_reduct@cell.embeddings)[2]
    colnames(sec_reduct@cell.embeddings) <- paste0(new_key, seq_len(ncols))
    # change assay used
    sec_reduct@assay.used <- assay_name
    # add reduction
    main_subset <- main_srt[, common_bcodes]
    main_subset[[new_name]] <- sec_reduct
    return(main_subset)
  }
  # involved case - we need to add NA for each cell in main_srt not present in
  # secondary_srt
  dropped_bcodes <- main_bcodes[which(!(main_bcodes %in% sec_bcodes))]
  sec_bcodes <- c(common_bcodes, dropped_bcodes)
  ncols <- dim(sec_reduct@cell.embeddings)[2]
  nrows <- length(sec_bcodes)
  # generate new cell embeddings
  new_embeddings <- matrix(data = rep(NA, each=nrows*ncols), nrow=nrows)
  colnames(new_embeddings) <- paste0(new_key, seq_len(ncols))
  rownames(new_embeddings) <- sec_bcodes
  new_embeddings[common_bcodes, ] <- as.matrix(
                                    sec_reduct@cell.embeddings[common_bcodes, ])
  # create dimensionality reduction object and add to Seurat
  new_reduct <- CreateDimReducObject(embeddings = new_embeddings,
                                     assay = assay_name,
                                     key = new_key)
  new_reduct@global <- sec_reduct@global
  main_srt[[new_name]] <- new_reduct
  return(main_srt)
}

#' Loads TCR clone information to a data.frame
#'
#' Inputs:
#' @param tcr_dir - the folder where TCR-seq data is stored. This folder should
#' contain both the "filtered_contig_annotations.csv" and "clonotypes.csv"
#' files from the same sample
#' @param fix_bcode - whether to remove the "-1" from 10x barcodes. Should
#' match what was used for loading GEx

#' Output:
#' @param clone_info - a data.frame with TCR information for each barcode a TCR
#' was detected in. Information includes everything in
#' filtered_contig_annotations, plus total counts in the experiment, a single
#' string version of the CDR3 AA sequences, and the clonality of the associated
#' clone
#'
load_vdj <- function(tcr_dir, fix_bcode=TRUE, verbose=FALSE) {
  # make sure this is loading a single dataset
  if (length(tcr_dir)!=1) {
    cat("Did not provide a single folder!\n")
    return(NULL)
  }
  if (!dir.exists(tcr_dir)) {
    cat(tcr_dir, "does not exist!\n")
    return(NULL)
  }
  tcr_file <- list.files(tcr_dir, pattern = "filtered.*csv",
                         full.names = T) #, recursive = T)
  if (length(tcr_file) == 0) {
    tcr_file <- list.files(paste(tcr_dir, "outs", sep = "/"),
                           pattern = "filtered.*csv",
                           full.names = T)
  }
  tcr_cells <- read.csv(tcr_file)

  # remove contigs that are not productive
  has_type <- grepl(pattern = "clonotype", x = tcr_cells$raw_clonotype_id)
  is_cell <- as.logical(tcr_cells$is_cell)
  confident <- as.logical(tcr_cells$high_confidence)
  is_fl <- as.logical(tcr_cells$full_length)
  is_productive <- as.logical(tcr_cells$productive)
  filter_tcrs <- has_type & is_cell & confident & is_fl & is_productive
  tcr_cells <- tcr_cells[filter_tcrs, ]

  if (fix_bcode) {
    tcr_cells$barcode <- gsub(pattern = "-1", replacement = "",
                              x = tcr_cells$barcode)
  }

  # add column for total counts; this will be helpful with merging experiments
  # and with calculating frequencioes
  total_counts <- length(unique(tcr_cells$barcode))
  tcr_cells[["total_counts"]] <- total_counts

  cdr3s_aa <- sapply(X = tcr_cells$raw_clonotype_id,
                     FUN = generate_cdr3s,
                     clone_info = tcr_cells)
  tcr_cells[["cdr3s_aa"]] <- cdr3s_aa

  clonalities <- sapply(X = tcr_cells$raw_clonotype_id,
                        FUN = calc_clonality,
                        clone_info = tcr_cells)
  tcr_cells[["TCR_clonality"]] <- clonalities

  if (verbose) {
    total_clones <- length(unique(tcr_cells$raw_clonotype_id))
    trb_rows <- tcr_cells$chain == "TRB"
    trbs <- unique(tcr_cells[trb_rows, "cdr3"])
    cat("Total clones detected:", total_clones,
        "\nUnique CDR3-B AA sequences:", length(trbs),
        "\n", sep=" ")
  }

  clone_info <- as.data.frame(tcr_cells)
  return(clone_info)
}


#' Add TCR information to a Seurat object. This function takes in a table of
#' the same format as filtered_contig_annotations and adds to metadata the
#' clonotype id of each barcode, the CDR3 AA sequences, and the clonality of
#' the corresponding TCR clone (either as number of cells, or fraction of total)
#'
#' Inputs:
#' @param srt_obj - a Seurat object to which to add the TCR information
#' @param clone_info - TCR information, in the format output by load_vdj
#' @param use_frequency - a boolean determining whether clonality metadata
#' should be the fraction of the population (TRUE) or absolute count (FALSE)
#'
#' Output:
#' @param srt_new - the same Seurat object with metadata added
#'
add_vdj <- function(srt_obj, clone_info,
                    use_frequency = TRUE, verbose=FALSE) {
  clonotypes = list("clonotype_id"=clone_info$raw_clonotype_id)
  clones_barcodes <- aggregate(x = clone_info$barcode,
                               by = clonotypes,
                               FUN = paste, collapse = ",")
  clonotype_ids <- clones_barcodes$clonotype_id
  clones_barcodes <- sapply(X = clones_barcodes[, 2],
                            FUN = function(bcodes) {
                              bcode_list <- strsplit(x=bcodes, split=",")[[1]]
                              bcode_list <- unique(bcode_list)
                              bcode_string <- paste(bcode_list, collapse = ",")
                              return(bcode_string)
                            }
                           )
  # get calculated clonality, total cells per group and CDR3 AA sequences. Done
  # by aggregation and taking first value, since they are unique per clonotype
  # and to avoid recalculating the same thing in ways that may not be identical
  columns <- c("cdr3s_aa", "TCR_clonality", "total_counts")
  clone_stats <- aggregate(x = clone_info[, columns],
                           by = clonotypes,
                           FUN = function(x) {x[1]})
  # generate Seurat metadata
  meta_data <- as.data.frame(x = list("TCR_clonality" = NULL,
                                      "clonotype_id" = NULL,
                                      "cdr3s_aa" = NULL))
  for (iclone in seq_len(length(clonotype_ids))) {
    clone_name <- clonotype_ids[iclone]
    barcodes_string <- clones_barcodes[iclone]
    clone_barcodes <- strsplit(x = barcodes_string,
                               split = ",")[[1]]
    nbarcodes <- length(clone_barcodes)
    clonality <- clone_stats$TCR_clonality[iclone]
    if (use_frequency) {
      # there may be multiple totals in data containing multiple samples
      total <- clone_stats$total_counts[iclone]
      clonality <- clonality / total
    }
    cdr3 <- clone_stats$cdr3s_aa[iclone]
    clone_meta <- list("TCR_clonality" = rep(x = clonality, each = nbarcodes),
                       "clonotype_id" = rep(x = clone_name, each = nbarcodes),
                       "cdr3s_aa" = rep(x = cdr3, each = nbarcodes))
    clone_meta <- as.data.frame(x = clone_meta, row.names = clone_barcodes)
    meta_data <- rbind(meta_data, clone_meta)
  }

  srt_new <- AddMetaData(object = srt_obj, metadata = meta_data)

  if (verbose) {
    total_clones <- length(clonotype_ids)
    trbs <- clone_info$cdr3[which(clone_info$chain == "TRB")]
    total_trb_num <- length(unique(trbs))

    gex_clones <- srt_new$cdr3s_aa
    gex_nclones <- length(unique(gex_clones))
    gex_trbs <- sapply(X = gex_clones, FUN = get_chain, chain = "TRB")
    gex_trbs <- unlist(sapply(X = gex_trbs,
                              FUN = function (x) {
                                unlist(strsplit(x=x, split=","))
                              }
                            ))
    gex_ntrbs <- length(unique(gex_trbs))

    cat("Of", total_clones, "total clones,",
        gex_nclones, "are retained in this data (",
        gex_ntrbs, "of", total_trb_num, "CDR3-B AA sequences)\n", sep=" ")
  }

  return(srt_new)
}


#' A helper function that generates the equivalent of the cdr3s_aa column of
#' clonotypes.csv from the entries of filtered_contig_annotations
#'
#' Inputs:
#' @param clone_info - a table or data.frame with the same entries as
#' filtered_contig_annotations/ the output of load_vdj
#' @param clone - the name of a clone, matching a subset of entries in
#' raw_clonotype_id
#'
#' Returns:
#' @param cdr3s_aa - a vector of strings representing all CDR3 AA sequences
#' detected for each clone. Each AA sequence is of the form "TR[AB]:####", where
#' ### is replaced by the AA sequence. Sequences from each detected chain are
#' separated by a semicolon
#'
generate_cdr3s <- function(clone, clone_info) {
  rows <- which(clone_info$raw_clonotype_id == clone)
  chains <- clone_info[rows, "chain"]
  cdr3s <- clone_info[rows, "cdr3"]
  cdr3_strs <- paste(chains, cdr3s, sep = ":")
  cdr3s_aa <- paste(sort(unique(cdr3_strs)), collapse = ";")
  return(cdr3s_aa)
}


#' A helper function that calculates the copy number of a given clone
#'
#' Inputs:
#' @param clone_info - a table or data.frame with the same entries as
#' filtered_contig_annotations/ the output of load_vdj
#' @param clone - the name of a clone, matching a subset of entries in
#' raw_clonotype_id
#'
#' Returns:
#' @param clonality - an integer representing the copy number of the clone
#'
calc_clonality <- function(clone, clone_info) {
  rows <- which(clone_info$raw_clonotype_id == clone)
  barcodes <- clone_info[rows, "barcode"]
  clonality <- length(unique(barcodes))
  return(clonality)
}


#' A function to merge 10X VDJ information from multiple samples. Typically,
#' these will be different sequencing runs from the same original cell
#' population. This function can deal with data with multiple subgroups of
#' runs. Within each such group, it will fully merge clones with the same
#' Beta-CDR3 AA sequence and join their clonality. It will not, however,
#' attempt to merge clones across defined groups.
#'
#' Inputs:
#' @param vdj_list - a list of objects containing TCR clone information. The
#' format of each entry should be the same as clone_info from load_vdj. The
#' list can be either named or unnamed.
#' @param sample_groups - a list, in which each entry is a list of the sample
#' numbers that are to be merged. If not specified, all samples will be merged.
#' Note that sample order should be the same order the samples were added to
#' the object in. Alternatively, these can be the names of the experiments
#' (in this case, these should be column names in vdj_list). If this list is
#' named, these names are treated as group names (relevant mostly if
#' merge_clone_id is TRUE).
#' @param merge_clone_id - a boolean argument. If TRUE, will rename all
#' clonotype IDs and merge those than have identical CDR3 AA sequences.
#' @param merge_clonality - a boolean argument. If TRUE, clonality is merged
#' (copy numbers of identical clones summed and total detected TCRs also summed)
#' across samples in the same group. Setting to FALSE means that, even if
#' identical clones get a common name, their clonality will still be in the
#' context of their individual samples.
#'
#' Output:
#' @param new_group_info - a combined object with the same format and combining
#' clone information for each sub-group. Barcodes are modified with the
#' corresponding column name, to make sure cell barcodes remain unique.
#'
merge_vdj <- function(vdj_list, sample_groups = NULL,
                      merge_clone_id = TRUE, merge_clonality = TRUE) {
  # sample names
  samples <- names(vdj_list)
  if (is.null(samples)) {
    samples <- seq_len(length(samples))
  }
  # manage sample_groups format
  if (is.null(sample_groups)) { # if not specified, merge all
    sample_groups <- list(seq.int(1, length(samples)))
  }
  # convert sample names to indexes
  if (typeof(sample_groups[[1]]) == "character") {
    temp_ <- vector(mode = "list", length = length(sample_groups))
    for (igr in seq_len(length(sample_groups))) {
      for (iexp in seq_len(length(sample_groups[[igr]]))) {
        temp_[[igr]][iexp] <- which(samples == sample_groups[[igr]][iexp])
      } #end iexp for
    } #end igr for
    # keep group names
    if (!is.null(names(sample_groups))) {
      names(temp_) <- names(sample_groups)
    }
    sample_groups <- temp_
    rm(temp_)
  }

  new_clone_info <- NULL
  # loop for merging TCR clones within each group
  for (isg in seq_len(length(sample_groups))) {
    new_group_info <- merge_vdj_group(vdj_list = vdj_list,
                                      sample_groups = sample_groups,
                                      samples = samples,
                                      isg = isg,
                                      merge_clone_id = merge_clone_id,
                                      merge_clonality = merge_clonality)
    new_clone_info <- rbind(new_clone_info, new_group_info)
  }
  return(new_clone_info)
}


#' Simple function that makes an empty data.frame with VDJ info. This returns
#' a data.fram with all filtered_contig_annotations entries, plus group-specific
#' information generated by merge_vdj
empty_vdj_group <- function() {
  new_clone_info <- as.data.frame(x =
                      list("barcode" = NULL, "is_cell" = NULL,
                           "contig_id" = NULL, "high_confidence" = NULL,
                           "length" = NULL, "chain" = NULL, "v_gene" = NULL,
                           "d_gene" = NULL, "j_gene" = NULL, "c_gene" = NULL,
                           "full_length" = NULL, "productive" = NULL,
                           "fwr1" = NULL, "fwr1_nt" = NULL, "cdr1" = NULL,
                           "cdr1_nt" = NULL, "fwr2" = NULL, "fwr2_nt" = NULL,
                           "cdr2" = NULL, "cdr2_nt" = NULL, "fwr3" = NULL,
                           "fwr3_nt" = NULL, "cdr3" = NULL, "cdr3_nt" = NULL,
                           "fwr4" = NULL, "fwr4_nt" = NULL, "reads" = NULL,
                           "umis" = NULL, "raw_clonotype_id" = NULL,
                           "raw_consensus_id" = NULL,
                           "total_counts" = NULL, "cdr3s_aa" = NULL,
                           "TCR_clonality" = NULL, "group" = NULL,
                           "sample" = NULL, "original_clone_id" = NULL))
  return(new_clone_info)
}

#' A list of columns to extract from 10X. Helps to keep things compatible
vdj_columns <- function() {
  c("barcode", "is_cell", "contig_id", "high_confidence", "length",
    "chain", "v_gene", "d_gene", "j_gene", "c_gene",
    "full_length", "productive",
    "cdr3", "cdr3_nt",
    "reads", "umis",
    "raw_clonotype_id", "raw_consensus_id",
    "total_counts", "cdr3s_aa", "TCR_clonality",
    "group", "sample", "original_clone_id")
}


#' A helper function to merge VDJ info within a group for which all information
#' is aggregated. See merge_vdj for more info.
merge_vdj_group <- function(vdj_list, sample_groups, samples, isg,
                            merge_clone_id = TRUE, merge_clonality = TRUE) {
  group_info <- NULL
  group_names <- names(sample_groups)
  if (is.null(group_names)) {
    group_name <- paste0("group", isg)
  }
  else {
    group_name <- group_names[isg]
  }

  # first concatenate all info
  total_count_list <- list()
  for (idx_ in sample_groups[[isg]]) {
    exp_name <- samples[idx_]
    exp_info <- vdj_list[[exp_name]]
    exp_info[["group"]] <- group_name
    exp_info[["sample"]] <- exp_name
    total_count_list[[exp_name]] <- exp_info[["total_counts"]][1]
    # modify barcodes and clonotypes to represent origin
    exp_info[["barcode"]] <- paste(exp_name, exp_info$barcode, sep="_")
    # for now, keep these. Will be updated if clones are merged
    exp_info[["raw_clonotype_id"]] <- paste(exp_name,
                                            exp_info$raw_clonotype_id,
                                            sep="_")
    exp_info[["original_clone_id"]] <- exp_info$raw_clonotype_id
    exp_info <- select(exp_info, vdj_columns())
    group_info <- rbind(group_info, exp_info)
  }

  # aggregate info by cdr3 AA sequence. Clones with multiple chains need to
  # have a complete match to aggregate.
  if (merge_clone_id) {
    aa_seqs <- list("aa_seqs" = group_info$cdr3s_aa)
    new_clones <- aggregate(x = group_info$raw_clonotype_id,
                            by = aa_seqs,
                            FUN = paste)
    new_cdr3s <- new_clones$aa_seqs
    old_clones <- new_clones[, 2]
    clone_counts <- seq_len(length(new_cdr3s))
    new_clone_names <- paste0(group_name, "_clonotype", clone_counts)
    # update raw_clonotype_id;
    # information on original clone name is kept in original_clone_id
    for (iclone in clone_counts) {
      rows <- which(group_info$raw_clonotype_id %in% old_clones[iclone][[1]])
      group_info[rows, "raw_clonotype_id"] <- new_clone_names[iclone]
    }
  }

  # merge clonalities; use clone name to allow simple separation of clonotype
  # and clonality merging
  if (merge_clonality) {
    group_info[["total_counts"]] <- sum(unlist(total_count_list))
    clonalities <- sapply(X = group_info$raw_clonotype_id,
                          FUN = calc_clonality,
                          clone_info = group_info)
    group_info[["TCR_clonality"]] <- clonalities
  }

  return(group_info)
}


#' Function to parse out chain (e,g, Beta-CDR3) sequences from the default
#' format of 10X VDJ data
get_chain <- function(cdr3, chain="TRB") {
    seqs <- strsplit(x = cdr3, split = ";")[[1]]
    chains <- grep(pattern = chain, x = seqs, value = T)
    max_len <- max(nchar(chains))
    chains <- substr(x = chains, start = 5, stop = max_len)
    chain_str <- paste(chains, collapse = ",")
    return(chain_str)
}


#' Generate a list of TCRs with a given characterisitc. Characteristics can
#' be defined by a set of possible value from a single metadata entry.
#'
#' Inputs:
#' @param srt_obj - a Seurat object from which to get the clones and
#' characteristics
#' @param target_meta - the characterisitc metadata name (e.g.
#' "seurat_clusters")
#' @param target_values - the possible values of target_meta.
#' Can be a single value or a vector
#' @param mode - can be eithr "select" or "range". In "select" mode,
#' target_meta values are expected to be present as a category in
#' target_values. In contrast, "range" mode interprets target_values as a
#' min-max value range, and target_values must have size 2.
#' @param get_meta - the metadata name for clonotype names. Can be modified
#' such that other type of metadata are aggregated besides TCR clones.
#' @param add_meta - a vector of additional metadata to be collected. Default
#' is to include clonality and CDR3 AA sequence.
#' @param funcs - functions for aggregating the data in add_meta. This needs
#' to be a list or vector of functions matching add_meta (e.g.
#' c(mean, unique), or function (x) {x[1]})
#' @param save_path - a path to which to save the resulting table, with
#' get_meta and add_meta from cells matching target_values. NULL means not to
#' save to disk (default)
#'
#' Returns:
#' @param final_targs - a data.frame with the aggregated data table
#'
gen_tcr_list <- function(srt_obj, target_meta, target_values, mode="select",
                         get_meta="clonotype_id",
                         add_meta=c("TCR_clonality", "cdr3s_aa"),
                         funcs=c(mean, unique),
                         save_path=NULL) {
  # identify target clones
  metadata <- srt_obj[[target_meta]][[1]]
  cells_clones <- srt_obj[[get_meta]][[1]]
  if (mode == "select") {  # discrete categories
    target_clones <- cells_clones[which(metadata %in% target_values)]
  } else if (mode == "range") {  # numerical range
    target_values <- force_range(target_values)
    if (is.null(target_values)) {
      return(NULL)
    }
    target_clones <- cells_clones[which(metadata >= target_values[1] &
                                        metadata <= target_values[2])]
  } else {  # illegal value
    cat("Illegal value passed for 'mode'. Choose 'select' or 'range'\n")
    return(NULL)
  }

  cells <- colnames(srt_obj)[which(cells_clones %in% target_clones &
                             !is.na(cells_clones))]
  target_table <- FetchData(srt_obj,
                            vars=c(get_meta, add_meta),
                            cells=cells)
  # aggregate data by clonotype
  clonotypes <- list()
  clonotypes[[get_meta]] <- target_table[[get_meta]]
  barcodes <- aggregate(x=rownames(target_table), by=clonotypes,
                        FUN=paste, collapse=",")
  final_targs <- list("barcodes"=barcodes[, 2])
  final_targs[[get_meta]] <- barcodes[, 1]
  for (imeta in seq_len(length(add_meta))) {
    meta <- add_meta[imeta]
    func <- funcs[[imeta]]
    data <- aggregate(x=target_table[[meta]], by=clonotypes, FUN=func)
    final_targs[[meta]] <- data[, 2]
  }
  final_targs <- as.data.frame(x=final_targs,
                               row.names=final_targs[[get_meta]])
  if (!is.null(save_path)) {
    write.csv(final_targs, save_path)
  }
  return(final_targs)
}

#' A helper function that makes sure a vector of length 2 (range) can be
#' interpreted as a mathematical range of numbers
#'
force_range <- function(range, default_min=0, default_max=Inf) {
  # check length
  if(length(range) != 2) {
    cat("Provided range must be of length 2!\n")
    return(NULL)
  }
  # make sure both entries are numeric
  if (is.null(range[1]) | is.na(range[1])) {
    range[1] <- default_min
  }
  if (is.null(range[2]) | is.na(range[2])) {
    range[2] <- default_max
  }
  # flip if necessary
  if (range[1] > range[2]) {
    range <- rev(range)
  }
  return(range)
}
