library(ggplot2)
args = commandArgs(T)

for ( each_file in args[1:length(args)]) {
	path_list = unlist(strsplit(each_file,"/")) #按照/拆分,
	file_name = path_list[length(path_list)]
	gene_name = unlist(strsplit(file_name,"\\."))[1] # 按照.拆分
	print(gene_name)
	#title_name = paste(gene_name, unlist(strsplit(file_name,"\\."))[2], sep=".")

	data = read.table(each_file, header=T)
	#data  = subset(data1,readsNumber>=3)
	p1 = ggplot(data,aes(x=readsNumber, y=..count..)) + geom_histogram(stat="count", color="black", fill="white") +
		xlab("supporting reads") + ylab("Molecule")+
		#scale_x_continuous(breaks=seq(0,1000,100), limits=c(0,1000)) + # 设置x轴的坐标范围（limits），和间隔（breaks）
		#scale_x_continuous(breaks=seq(0,100,10), limits=c(0,100)) + # 设置x轴的坐标范围（limits），和间隔（breaks）
		ggtitle(gene_name) + theme(plot.title = element_text(hjust = 0.5))  # 设置标题，并居中放置在图案上方。

	p2 = ggplot(data,aes(x=readsNumber, y=..count..)) + geom_histogram(stat="count", color="black", fill="white") +
		xlab("supporting reads") + ylab("Molecule")+
		scale_x_continuous(breaks=seq(0,1000,100), limits=c(0,1000)) + # 设置x轴的坐标范围（limits），和间隔（breaks）
		#scale_x_continuous(breaks=seq(0,100,10), limits=c(0,100)) + # 设置x轴的坐标范围（limits），和间隔（breaks）
		ggtitle(gene_name) + theme(plot.title = element_text(hjust = 0.5)) # 设置标题，并居中放置在图案上方。
	
	p3 = ggplot(data,aes(x=readsNumber, y=..count..)) + geom_histogram(stat="count", color="black", fill="white") +
		xlab("supporting reads") + ylab("Molecule")+
		#scale_x_continuous(breaks=seq(0,1000,100), limits=c(0,1000)) + # 设置x轴的坐标范围（limits），和间隔（breaks）
		scale_x_continuous(breaks=seq(0,100,10), limits=c(0,100)) + # 设置x轴的坐标范围（limits），和间隔（breaks）
		ggtitle(gene_name) + theme(plot.title = element_text(hjust = 0.5))  # 设置标题，并居中放置在图案上方。
	
	ggsave(paste(gene_name,"umicount.histogram.raw.png",sep="."),p1)
	ggsave(paste(gene_name,"umicount.histogram.filter.png",sep="."),p2)
	ggsave(paste(gene_name,"umicount.histogram.zoom.png",sep="."),p3)
}
