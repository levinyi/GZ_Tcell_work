#!/usr/bin/env Rscript

# This is a script meant to automate analysis of 10X output, generate
# some basic plots and save output files and figures for further, more detailed
# analysis (e.g. using BBrowser). Assumes data has been through CellRanger.
# All downstream analysis is performed using the R package Seurat.
#
# Author: Seon Kinrot (RootPath Genomics)

# imports and requirements
suppressMessages(library(here))
suppressMessages(library(optparse))
suppressMessages(library(rBCS))
suppressMessages(source(here("src", "AutomationFunctions.R")))
suppressMessages(source(here('src', 'SeuratWrappers.R')))

# define CLI interface and parser
nl <- "\n\t\t"
parser <- OptionParser(prog="BatchSeurat",
                       description=paste("", "Command line interface for",
                                         "automated batch Seurat analysis"),
                       usage=paste("usage: %prog [-h] [-g] [-c | -s] [-m] [-v]",
                                   "[-i INPUT] [-o OUTPUT]",
                                   "[sample 1] [sample 2]...", sep = "\n"),
                       epilogue=paste0("For support, contact Seon Kinrot ",
                                       "(RootPath Genomics)"),
                       formatter=IndentedHelpFormatter)
parser <- add_option(parser, c("-g", "--gex_only"),
                     action="store_false", dest="has_vdj", default=TRUE,
                     help=paste("Only analyze GEX data, ignoring VDJ.",
                                "Note this has no bearing on CSP flags",
                                sep = nl))
parser <- add_option(parser, c("-c", "--cite_seq"),
                     action="store_true", dest="has_csp", default=FALSE,
                     help=paste("Specify this data includes CSP.",
                                "By default, assumes CSP and GEX are fused",
                                "(i.e. 'sep_csp' is False)", sep = nl))
parser <- add_option(parser, c("-s", "--sep_csp"),
                     action="store_true", default=FALSE,
                     help=paste0("Specify that GEX and CSP matrices are ",
                                 "separate.", nl, "If this flag is passed, ",
                                 "'cite_seq' is set to True"))
parser <- add_option(parser, c("-m", "--mouse"),
                     action="store_true", default=FALSE,
                     help=paste0("Specify this is a mouse dataset ",
                                 "[default is human]"))
parser <- add_option(parser, c("-v", "--verbose"),
                     action="store_true", default=FALSE,
                     help=paste("Print status updates as analysis advances",
                                "[default is not to print updates]",
                                sep = nl))
parser <- add_option(parser, c("-i", "--input"), default=getwd(),
                     help=paste0("Set the path to the input directory ",
                                 "containing the 10X analysis", nl,
                                 "[default is the current working ",
                                 "directory]"))
parser <- add_option(parser, c("-o", "--output"),
                     help=paste("Set the path to which output files are",
                                "saved [default is 'input'/analysis]",
                                sep = nl))

parsed_args <- parse_args2(parser)
in_list <- handle_options(parsed_args$options,
                          args = c("has_vdj", "has_csp", "sep_csp",
                                   "mouse", "verbose", "input", "output"))
in_list <- handle_pos_args(parsed_args$args, in_list)

# error checks
if (is.null(in_list)) {
  # error message should have already been printed from handle_cli
  quit(save="no", status=1, runLast=FALSE)
}

verbose <- in_list$verbose
if (verbose) {
  cat("Starting analysis!\n")
  start_time <- Sys.time()
}

# load data
data_list <- load_all_data(in_list)
if (verbose) {
  cat("Loaded raw data.",
      paste0("Elapsed time: ", time_string(start_time)),
      "Performing QC...\n",
      sep = "\n")
}

# Interactive/ Shiny QC
shiny_dir <- here("ShinyQC")
data_list$gex_list <- run_qc_shiny(host="192.168.2.185",data_list$gex_list,
                                   shiny_dir,
                                   in_list,
                                   port=8886)
if (verbose) {
    qc_time <- Sys.time()
    cat("QC completed. Integrating data...\n")
}

# data integration
n_projects <- length(in_list$project_names)
# edge case - no integration needed
if (n_projects == 1) {
   integrated_gex <- data_list$gex_list[[1]]
   integrated_csp <- data_list$csp_list[[1]]
   merged_vdj <- data_list$vdj_list[[1]]
} else {
   # integrate GEX
   integrated_gex <- suppressWarnings(
                        integrate_datasets(
                          srt_list = data_list$gex_list,
                          normalization_method = "SCT",
                          regress = c("pct_mito", "nCount_RNA"),
                          nfeatures = 3000,
                          rmv_grep = "TR[AB]V",
                          readd_counts = TRUE
                        )
                     )

   # integrate CSP
   integrated_csp <- NULL
   if (in_list$has_csp) {
     integrated_csp <- suppressWarnings(
                          integrate_datasets(
                            srt_list = data_list$csp_list,
                            normalization_method = "CLR",
                            nfeatures = 2000,
                            readd_counts = TRUE
                          )
                        )
   }

   # VDJ merging
   merged_vdj <- NULL
   if (in_list$has_vdj) {
     # Only supports merging of all samples
     merged_vdj <- suppressWarnings(
                      merge_vdj(vdj_list = data_list$vdj_list)
                   )
   }
}
if (verbose) {
   cat("Data integration complete.",
       paste0("Integration runtime: ", time_string(qc_time)),
       paste0("Total runtime: ", time_string(start_time)),
       "Clustering and plotting...\n",
       sep = "\n")
}

 # cluster independently and combine all data into GEX Seurat object
if (in_list$has_csp) {
   # add to GEX Seurat object; this also clusters CSP
   integrated_gex <- suppressWarnings(
                        fuse_cite_seq(integrated_gex,
                                      integrated_csp,
                                      integrated = n_projects > 1
                                    )
                     )
}
# cluster GEX, but limit to cells that have CSP
if (n_projects == 1) {
    integrated_gex <- suppressWarnings(
                        normalize_and_find_var_feats(
                          integrated_gex,
                          normalization_method = "SCT",
                          regress = c("pct_mito", "nCount_RNA"),
                          rmv_grep = "TR[AB]V"
                        )
                      )
}
integrated_gex <- suppressWarnings(cluster_seurat(srt_obj = integrated_gex))
# add TCRs as metadata
if (in_list$has_vdj) {
    integrated_gex <- suppressWarnings(
                        add_vdj(srt_obj = integrated_gex,
                                clone_info = merged_vdj,
                                use_frequency = TRUE
                              )
                      )
}
# generate plots
integrated_gex <- suppressMessages(batch_plots(integrated_gex, in_list))
if (verbose) {
    cat("Analysis complete!",
        paste0("Total runtime: ", time_string(start_time)),
        "Saving data to disk...\n",
        sep = "\n")
}

# save analyzed data
if (in_list$has_vdj) {
    write.csv(merged_vdj,
              file = paste0(in_list$save_dir, "/merged_tcr_table.csv")
            )
    # generate table to assist TCR selection
    feature_list <- NULL
    assay_list <- NULL
    if (in_list$has_csp) {
      csp_features <- c("PD1-totSeq", "CD39-totSeq", "CD69-totSeq",
                        "TCRAV7.2-totSeq", "CD25-totSeq", "CD4-totSeq",
                        "CD8A-totSeq", "ITGAE-totSeq")
      feature_list <- c(feature_list, csp_features)
      csp_assays <- rep("ADT", length.out = length(csp_features))
      assay_list <- c(assay_list, csp_assays)
    }
    metadata_list <- paste0(names(GLOBAL_MODULES), "_AUCscore")
    save_path <- paste0(in_list$save_dir, "/tcr_stats_barcodes.csv")
    selection_table <- suppressMessages(
                          tcr_selection_table(vdj_data = merged_vdj,
                                              srt_obj = integrated_gex,
                                              features = feature_list,
                                              assays = assay_list,
                                              metadata = metadata_list,
                                              save_path = save_path
                                            )
                        )
}
if (in_list$has_csp) {
    saveRDS(integrated_csp,
            file = paste0(in_list$save_dir, "/integrated_csp.rds"))
}
saveRDS(integrated_gex,
        file = paste0(in_list$save_dir, "/integrated_seurat.rds"))
# export to bcs format for BBrowser
gex_assay <- "integrated"
csp_assay <- "adt_integrated"
if (n_projects == 1) {
    gex_assay <- "SCT"
    csp_assay <- "ADT"
}
a <- suppressMessages(
    rBCS::ExportSeurat(integrated_gex,
                       paste0(in_list$save_dir, "/integrated_data_BBrowser.bcs"),
                       author = "SK",
                       raw.rna = gex_assay,
                       norm.rna = gex_assay,
                       raw.adt = csp_assay,
                       norm.adt = csp_assay)
)
if (verbose) {
    cat("All done! Have a great day!\n")
}
