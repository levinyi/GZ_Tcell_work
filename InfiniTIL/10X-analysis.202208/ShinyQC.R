#!/usr/bin/env Rscript

# This is a wrapper script to run the Shiny App ShinyQC, which generates
# a UI for visual QC of 10X RNA sequencing data
#
# To run this script, make sure to save input data (a list of Seurat objects)
# in the ShinyQC/data folder, under the name "input.rds". After QC filters are
# chosen and the App is closed, the filtered Seurat objects will be saved as
# a list in the same order as "output.rds" in the same directory.

# Author: Seon Kinrot (RootPath Genomics)

suppressMessages(library(here))
suppressMessages(library(shiny))
suppressMessages(library(optparse))

parser <- OptionParser(prog="ShinyQC",
                       description=paste("", "Command line interface for",
                                         "semi-automated 10X QC filtering"),
                       usage="usage: %prog [-h] [-p PORT]",
                       epilogue=paste0("For support, contact Seon Kinrot ",
                                       "(RootPath Genomics)"),
                       formatter=IndentedHelpFormatter)
parser <- add_option(parser, c("-p", "--port"), default=8886,
                    help=paste("Set the port to use for the Shiny App",
                               "(default: 8886)"))
args <- parse_args(parser)

shiny_dir <- here("ShinyQC")
runApp(shiny_dir, launch.browser = TRUE, port = args$port)
stopApp() # shouldn't be necessary, but just in case
