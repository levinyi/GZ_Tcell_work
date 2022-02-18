# data1 = "/cygene2/work/P0000-Blackbird/2103-GA001/GA001003/analysis/GA001003_WES_RNAseq/GA001003.variants.funcotated.with.minigene.MAF.add.TPM.RNA_Depth.xls"
# data2 = "GA13_TU_CL1-P3-gDNA.variants.funcotated.with.minigene.MAF.xls"

library(dplyr)
library(ggplot2)
args = commandArgs(T)
#args = ["/cygene2/work/P0000-Blackbird/2103-GA001/GA001003/analysis/GA001003_WES_RNAseq/GA001003.variants.funcotated.with.minigene.MAF.add.TPM.RNA_Depth.xls",
#	"GA13_TU_CL1-P3-gDNA.variants.funcotated.with.minigene.MAF.xls",
# 	"GA13", 
# 	"GA13_TU_CL1-P3"]

data1 = args[1]
data2 = args[2]
suffix_left = args[3]
suffix_right = args[4]

table1 = read.table(data1, header=T)
table2 = read.table(data2, header=T)

full_table <- full_join(table1, table2, by=c("Hugo_Symbol","Start_Position"),suffix = c(paste('.',suffix_left,sep=''),paste('.',suffix_right,sep='')))
write.table(full_table, file=paste(suffix_left,"_VS_", suffix_right,".Raw.xls", sep=""), sep="\t", row.names = FALSE, quote=FALSE)

selected_table <- full_table %>% select(Hugo_Symbol,paste('tumor_f.',suffix_left,sep=''),paste('tumor_f.',suffix_right, sep='')) %>% na.omit()
write.table(selected_table, file=paste(suffix_left,"_VS_", suffix_right,".subset.VAF.xls",sep=""), sep="\t", row.names = FALSE, quote=FALSE)

# correlation diagram
# ggplot(select_table) + geom_point() 

# venn diagram

