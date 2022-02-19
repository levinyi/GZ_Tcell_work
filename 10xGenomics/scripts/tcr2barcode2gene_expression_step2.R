library(Seurat)
args= commandArgs(T)

RNAseq_path_10x = args[1]
# e.g. P0000-blackbird/2103-HC001-ZJW/G85E2L1_RNA-Seq/filtered_feature_bc_matrix
candidate_barcode = args[2] 
# e.g. "total.TCR.clonotype.in.raw_barcode.csv"
outputfile = args[3]

HC = Read10X(RNAseq_path_10x)
HC = CreateSeuratObject(counts = HC, min.cells = 3)
###################
# read TCR data
TCR_barcode_info = read.csv(candidate_barcode)
HC@meta.data = cbind(HC@meta.data, TCR_barcode_info[,-1])

#######################################################################################

###########################################################################################
cell.use = rownames(HC@meta.data[which(HC@meta.data$expansion =="Known"),])
sub_HC = subset(HC, cells =cell.use)


# convert a sparse matrix to a data frame like this:
gene_expression_matrix = as.data.frame(sub_HC@assays$RNA@counts)
# filter gene by 0 line.
# gene_expression_matrix_filtered = gene_expression_matrix[which(rowSums(gene_expression_matrix)>0),]
gene_expression_matrix_filtered = gene_expression_matrix
# rbind barcode info to matrix. 
info.use = TCR_barcode_info[which(TCR_barcode_info$expansion =="Known"),]
info.use.t = t(info.use)
colnames(info.use.t) =  info.use.t[1,]

rbind_info = rbind(info.use.t, gene_expression_matrix_filtered, stringsAsFactors = FALSE)

# # useless
library(tibble)
write.table(rbind_info[-c(1,4),] %>% rownames_to_column("barcode"), file = outputfile, sep = "\t", row.names =FALSE, quote = FALSE)
