#' This script is used after scRNAseq Infinitil pipeline.
setwd("/cygene/work/00.test/pipeline/10xGenomics/scripts")

# generate a table

library(Seurat)

#' input: tcr_stats_barcodes.csv
#' output: table, unique clonotypes with cell-subtype information from GXP of CBs and FACS annotation.
#' 
#' 
FOXP3 = 20

tcr_stats = read.csv("/cygene2/pipeline/10X/data/SA501001-TIL02/analysis/tcr_stats_barcodes.csv")
tcr_stats[1:4,]
tcr_stats %>% select(c("renamed_clone","TCR_clonality"))
