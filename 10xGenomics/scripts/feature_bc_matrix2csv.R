#!/usr/bin/env Rscript

# USAGE:
# Rscript feature_bc_matrix2csv.R MATRIX_PATH
# example:
# Rscript feature_bc_matrix2csv.R G55E2L1/outs/filtered_feature_bc_matrix

package_list <- c("Seurat")
for(p in package_list){
  if(!suppressWarnings(suppressMessages(require(p, character.only = TRUE, quietly = TRUE, warn.conflicts = FALSE)))){
    install.packages(p, repos="http://cran.r-project.org")
    suppressWarnings(suppressMessages(library(p, character.only = TRUE, quietly = TRUE, warn.conflicts = FALSE)))
  }
}

#####################################################

library(Seurat)
args = commandArgs(T)

RNA_seq = Read10X(args[1])

genes = c("PDCD1","CD4","CD8A","CD8B")

RNA_seq_select = RNA_seq[(rownames(RNA_seq) %in% genes ),]
RNA_seq_select_T = t(as.matrix(RNA_seq_select))
write.table(RNA_seq_select_T, file=paste(args[1],"Gene_expressed.count.csv",sep="/"),sep=",", row.names = TRUE, col.names=NA, quote=FALSE)














