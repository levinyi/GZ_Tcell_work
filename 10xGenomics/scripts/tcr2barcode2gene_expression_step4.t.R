args= commandArgs(T)
input = args[1]
outputfile = args[2]

# input = "/cygene2/work/P0000-Blackbird/2103-BL001/BL001001/BL001001_G178E2L1_scRNAseq/outs/candidate.tcr.gene_expression.matrix.raw.marked.subset.xls"

f = read.csv(input, header=F,sep = "\t")
f = t(f)
# outputfile = "/cygene2/work/P0000-Blackbird/2103-BL001/BL001001/BL001001_G178E2L1_scRNAseq/outs/test.xls"

write.table(f, file = outputfile, sep = "\t", row.names =FALSE, col.names = FALSE, quote = FALSE)
