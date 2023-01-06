#' This script is used after scRNAseq Infinitil pipeline.


# generate a table

library(Seurat)

#' input: tcr_stats_barcodes.csv
#' output: table, unique clonotypes with cell-subtype information from GXP of CBs and FACS annotation.
#' 
#' 

tcr_stats = read.csv("tcr_stats_barcodes.csv")
