library(tidyr)
library(dplyr)
library(ggplot2)

# args = [GA13_TU_CL1-P3-gDNA.Tumor.bam.DepthOfCoverage.txt]
# this script is only suit for GATK DepthOfCoverage result.

args = commandArgs(T)

for (each_file in args){
	output = dirname(each_file)
	file_name = basename(each_file)
	title_name = paste(unlist(strsplit(file_name, "\\."))[1:2], collapse=".")
	print(title_name)
	
	data = read.csv(each_file, header=T)

	data2 = data %>% separate(Locus,c("chr","pos"),":") %>%
		# separate(chr,c("chr1","chr"),"chr") %>%  
		group_by(chr) %>%
		summarise(chr_len = n(), 
			  total_depth = sum(Total_Depth),
			  coverage_leng = sum(Total_Depth!=0)) %>%
		mutate(rate1 = total_depth/chr_len, rate2=coverage_leng/chr_len)
	write.table(data2, file=paste(paste(output, title_name,sep = "/"),"DepthOfCoverage.rate.xls",sep="."), 
	            sep="\t", row.names = FALSE, quote=FALSE)
}
print("all done")
