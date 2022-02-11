library(tidyr)
library(dplyr)
library(ggplot2)

# args = [GA13_TU_CL1-P3-gDNA.Tumor.bam.DepthOfCoverage.txt]
args = commandArgs(T)


for (each_file in args){
	path_list = unlist(strsplit(each_file,"/"))
	file_name = path_list[length(path_list)]
	title_name = paste(unlist(strsplit(file_name, "\\."))[1:2], collapse=".")
	print(title_name)
	data = read.csv(each_file, header=T)

	data3 = data %>% separate(Locus,c("chr","pos"),":") %>%
		# separate(chr,c("chr1","chr"),"chr") %>%  
		group_by(chr) %>%
		summarise(chr_len = n(), 
			  total_depth = sum(Total_Depth),
			  coverage_leng = sum(Total_Depth!=0)) %>%
		mutate(rate1 = total_depth/chr_len, rate2=coverage_leng/chr_len)
	write.table(data3, file=paste(title_name,"coverage.depth.rate.xls",sep="."), sep="\t", row.names = FALSE, quote=FALSE)
}
