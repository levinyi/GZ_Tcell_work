library(ggplot2)
args = commandArgs(T)
data = read.table(args[1], header = F)
name = args[2]
ggplot(data, aes(x=reorder(data$V1,-data$V4), y=reorder(data$V2,data$V5), fill=log(data$V3,2))) + 
	geom_raster() + ggtitle(name)+
	xlab("UMIs") + ylab("clonotype")+
	theme(legend.title = element_blank(),
	      axis.text = element_blank(),
	      )
output = paste(name,".umi.raster.jpg",sep=".")
ggsave(output)
print(output)
