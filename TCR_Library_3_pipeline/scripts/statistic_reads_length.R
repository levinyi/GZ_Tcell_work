library(ggplot2)
args = commandArgs(T)
result_dir = args[1]
for (each_file in args[2:length(args)]){
	a = unlist(strsplit(each_file,"/"))
	name = a[length(a)]
	name = unlist(strsplit(name,"\\."))[1]
	print(name)
	data = read.table(each_file, header=F)

	ggplot(data,aes(x=V2,y=..count..)) + 
		geom_histogram(stat="bin",binwidth = 1,color="black", fill="white")+
		scale_x_continuous(breaks=seq(0,500,50),limits=c(0,500))+
		xlab("reads length")
	output_name = paste(name,"reads.cuted.length.jpg",sep='.')
	ggsave(paste(result_dir,output_name,sep="/"))
	print( paste(name,"reads.cuted.length.jpg",sep="."))

}
