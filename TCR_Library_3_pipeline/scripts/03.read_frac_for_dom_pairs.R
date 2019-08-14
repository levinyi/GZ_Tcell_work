library(ggplot2)
require(cowplot)

args = commandArgs(T)
result_dir = args[1]
for (each_file in args[2:length(args)]) {
	a =  unlist(strsplit(each_file, "/"))
	name = a[length(a)]
        name = unlist(strsplit(name, "\\."))[1]

        data = read.table(each_file, header=F)
	data2 = data[data$V3=='frac',]
	data3 = data[data$V3=='total',]

	p1 = ggplot(data3,aes(x=V4,y=V5,fill=V3,color=V2)) + geom_line() +
		scale_x_continuous(limits = c(0,1),breaks=seq(0,1,0.2)) +
		xlab("Thredhold for dominance") +ylab("# of Dominant pairs")+ 
		ggtitle(name) + theme(plot.title = element_text(hjust = 0.5)) +
		theme(legend.title = element_blank())

	p2 = ggplot(data2, aes(x=V4,y=V5,fill=V3,color=V2)) + geom_line() +
		scale_x_continuous(limits = c(0,1),breaks = seq(0,1,0.2)) +
		scale_y_continuous(limits = c(0,1),breaks = seq(0,1,0.2)) +
		xlab("Thredhold for dominance") + ylab("Frac. reads from dom. pairs") +
		ggtitle(name) + theme(plot.title = element_text(hjust = 0.5)) +
		theme(legend.title = element_blank()) # remove legend title
    	ggdraw() + 
		draw_plot(p1, 0,.5, 1, .5) +
		draw_plot(p2, 0, 0, 1, .5) +
		draw_plot_label(c("A","B"),c(0,0),c(1,.5),size=15)
	output_name = paste(name, 'dominant_pairs_graph.jpg',sep='.')
        ggsave(paste(result_dir,output_name,sep="/"),width=10,height=10)
	print("finished read_frac_for_dom_pairs.R")


}

