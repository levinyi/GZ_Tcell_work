library(Seurat)
args = commandArgs(T)
# RNA_seq = Read10X("G55E2L1/outs/filtered_feature_bc_matrix")
RNA_seq = Read10X(args[1])
#print(RNA_seq[1:10,1:6])

genes = c("PDCD1","CD4","CD8A","CD8B")

RNA_seq_select = RNA_seq[(rownames(RNA_seq) %in% genes ),]
# RNA_seq_select = RNA_seq[(rownames(RNA_seq) %in% genes ),1:9]
RNA_seq_select_T = t(as.matrix(RNA_seq_select))
write.table(RNA_seq_select_T, file=paste(args[1],"Gene_expressed.count.csv",sep="/"),sep=",", row.names = TRUE, col.names=NA, quote=FALSE)














