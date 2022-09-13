#!/usr/bin/env Rscript

# This is a wrapper script to run the CellRanger portion of BatchAnalysis
# without downstream Seurat analysis.
#
# To run this script, make sure all fastq files are saved in the same folder,
# and follow a predictable naming pattern. Names should include the sample type
# (GEX, VDJ or CSP), and files related to the same sample should share
# identical names (aside from sample type).
#
# Author: Seon Kinrot (RootPath Genomics)

# imports and requirements
suppressMessages(library(here))
suppressMessages(library(optparse))
suppressMessages(source(here("src", "AutomationFunctions.R")))

# define CLI interface and parser
nl <- "\n\t\t"
parser <- OptionParser(prog="BatchCR",
                       description=paste("", "Command line interface for",
                                         "automated CellRanger count runs"),
                       usage=paste("usage: %prog [-h] [-u CPUS] [-x MAX_MEM]",
                                   "[-g] [-c | -s] [-m] [-v]",
                                   "-t TRANSCRIPTOME [-j VDJ_REF]",
                                   "[-f FEATURE_REF] [-i INPUT]",
                                   sep = "\n"),
                       epilogue=paste0("For support, contact Seon Kinrot ",
                                       "(RootPath Genomics)"),
                       formatter=IndentedHelpFormatter)
parser <- add_option(parser, c("-u", "--cpus"), default=24,
                     help=paste("The number of cores to use for CellRanger",
                                "(default 24)"))
parser <- add_option(parser, c("-x", "--max_memory"), default=256,
                     help=paste("Max memory usage in CellRanger, in GB",
                                "(default 256"))
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
parser <- add_option(parser, c("-t", "--transcriptome"),
                     help=paste("The path to the transcriptome reference",
                                "file. *This is a required argument*",
                                sep = nl))
parser <- add_option(parser, c("-j", "--vdj_ref"),
                     help=paste("The path to the VDJ reference file",
                                "(required unless the '-g' flag is passed",
                                sep = nl))
parser <- add_option(parser, c("-f", "--feature_ref"),
                     help=paste("The path to the feature reference file",
                                "(required if passing the '-c' or '-s' flags)",
                                sep = nl))
parser <- add_option(parser, c("-i", "--input"), default=getwd(),
                     help=paste0("Set the path to the input directory ",
                                 "containing the 10X analysis", nl,
                                 "[default is the current working ",
                                 "directory]"))
# this is necessary to conform to handle_cli
parsed_args <- parse_args(parser)

# parse input
in_list <- handle_options(parsed_args,
                          args = c("cpus", "max_memory", "has_vdj",
                                   "has_csp", "sep_csp", "mouse", "verbose",
                                   "transcriptome", "vdj_ref", "feature_ref",
                                   "input"))
in_list <- handle_pos_args(list(), in_list)

# error checks
if (is.null(in_list)) {
# error message should have already been printed from handle_cli
quit(save="no", status=1, runLast=FALSE)
}
# check we have everything we need to run CellRanger
validated <- validate_cr(in_list)
if (!validated) {
  # error message should be printed from validate_cr
  quit(save="no", status=1, runLast=FALSE)
}

# time to get going!
verbose <- in_list$verbose
setwd(in_list$data_dir)
if (verbose) {
  cat("Running CellRanger...\n")
  cr_start <- Sys.time()
}

# start loop
n_projects <- length(in_list$project_names)
for (iproject in seq_len(n_projects)) {
  if (verbose) {
    cat("Running sample", in_list$project_names[[iproject]], "\n", sep = " ")
    sample_start <- Sys.time()
  }

  invisible(run_cr_sample(in_list, iproject))

  if (verbose) {
    cat(paste0("Done with sample ", in_list$project_names[[iproject]]),
        paste0("Time for this sample: ", time_string(sample_start), "\n"),
        sep = "\n")
  }
}  # end iproject loop
if (verbose) {
  cat("Done with CellRanger runs!",
      paste0("Total tiume: ", time_string(cr_start), "\n"),
      sep = "\n")
}
