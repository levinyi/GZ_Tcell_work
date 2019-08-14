library(ggplot2)
args = commandArgs(T)

result_dir = args[1]
for (each_file in args[2:length(args)]){
	a = unlist(strsplit(each_file,"/"))
	name = unlist(strsplit(a[length(a)],"\\."))[1]
	d = read.table(each_file,header=T)
	p = ggplot(d,aes(x=d$Rank,y=d$Value,color=d$CAT)) + 
		geom_point(shape=1,size=log(d$Reads)) +theme(panel.grid =element_blank()) +
		theme(axis.line = element_line(size=1, colour = "black")) + # # 再加上坐标轴（无刻度、无标签）
		xlab("ClonoType alpha") + ylab("Fraction Dominant") +
		ggtitle(name)+
		theme(plot.title = element_text(hjust = 0.5),legend.title=element_blank())

	output_name = paste(name,"top_pair_dominant.jpg",sep=".")
	ggsave(paste(result_dir,output_name,sep="/"))
	print(paste(output_name,"was finished draw!",sep=" "))
}
