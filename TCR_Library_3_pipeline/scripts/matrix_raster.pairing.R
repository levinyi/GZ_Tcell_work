if (!suppressWarnings(require("ggplot2"))){
 install.packages("ggplot2")
 require("ggplot2")
}
#library(dplyr)

args = commandArgs(T)
result_dir = args[1]
for (each_file in args[2:length(args)]) {
	path_list = unlist(strsplit(each_file,"/"))
	file_name = path_list[length(path_list)]
	title_name = unlist(strsplit(file_name, "\\."))[1]

	data = read.table(each_file, header=F)
	names(data) = c("a","b","freq")
	#print(head(data))
	data1 = head(data[order(data$freq,decreasing=T),],500)
	#data1 = arrange(data,desc(data$freq))
	p = ggplot(data1, aes(x=reorder(a,-freq),y=reorder(b,freq),fill=freq)) + geom_raster() + 
	    scale_fill_gradient(low="blue", high="red") + 
#	    scale_fill_gradientn(colours=c("#0000FFFF","#FFFFFFFF","#FF0000FF")) +
	    theme(axis.text=element_blank(),axis.ticks=element_blank(),legend.title=element_blank()) + 
	    ggtitle(title_name)+ 
	    theme(plot.title = element_text(hjust = 0.5)) +
	    xlab("TRA clones") + ylab("TRB clones")
   	outputname = paste(title_name, "pairing.matrix.jpg", sep = ".")
	ggsave(paste(result_dir,outputname,sep="/"),width=10,height=10)
	print("finished matrix_raster.pairing.R!")
}
