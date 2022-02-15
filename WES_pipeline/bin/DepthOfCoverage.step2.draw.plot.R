library(tidyr)
library(dplyr)
library(ggplot2)

# args = [GA13_TU_CL1-P3-gDNA.Tumor.coverage.depth.rate.xls]
# chr     chr_len total_depth     coverage_leng   rate1   rate2
# chr1    4046489 234630948       3795965 57.983834380867 0.938088550345744
# chr10   1488150 86518107        1464131 58.138028424554 0.983859825958405

args = commandArgs(T)


theme_x_axis = theme(
  axis.text.x = element_text(size = 15),
  axis.text.y = element_text(size = 15),
)

common_themes = theme(plot.title = element_text(size = 17, hjust = 0.5),
                      panel.grid.major = element_line(colour = NA), # 去掉网格线
                      panel.background = element_blank(), # 去掉背景
                      panel.grid.major.x = element_line(),  # 横向网格线linetype=2,xuxian
                      panel.grid.major.y = element_line(), # 纵向网格线
                      panel.border = element_rect(fill=NA), # 边框
                      axis.line = element_line(), # 坐标轴线
                      axis.title = element_text(size = 17),
                      legend.title = element_blank(),
                      legend.position = c(0.95,0.95), # bottom, right, 'left'
                      legend.background = element_rect(colour = NA, fill = NA),
)

for (each_file in args){
	path_list = unlist(strsplit(each_file,"/"))
	file_name = path_list[length(path_list)]
	title_name = paste(unlist(strsplit(file_name, "\\."))[1:2], collapse=".")
	print(title_name)
	data = read.table(each_file, header=T)
	# print(data)
	p = ggplot(data) + geom_line( aes(x=chr, y=rate2*100/1),group=1) +
		geom_point(aes(x=chr, y=rate2*100/1)) +
		geom_bar(aes(x=chr, y=rate1), stat = 'identity',fill = "#8CC53E", width = 0.7) +
		xlab("")+
		scale_y_continuous(name=expression("Mean Depth"), 
				   limits = c(0, 120),
				   sec.axis = sec_axis(~./100, name = "Proportion of Covered Bases")) +
		theme_bw() +theme_x_axis + common_themes +
		theme(axis.text.x = element_text(angle = 90, hjust=0.1, vjust=0.5))
	ggsave(filename = paste(title_name, "coverage.png", sep="."), plot = p, width=9, height=7, path = "./" )
	#write.table(data3, file=paste(title_name,"coverage.depth.rate.xls",sep="."), sep="\t", row.names = FALSE, quote=FALSE)
	print(paste(title_name,"coverage.png", sep="."))
}
