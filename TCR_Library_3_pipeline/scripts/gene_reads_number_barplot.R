library(ggplot2)
require(cowplot)
library(RColorBrewer)

args = commandArgs(T)
plot_dir = args[1]

for (each_file in args[2:length(args)]) {
    path_list = unlist(strsplit(each_file,"/"))
    file_name = path_list[length(path_list)]
    sample_name = unlist(strsplit(file_name, "\\."))[1]
    type_name = unlist(strsplit(file_name, "\\."))[2]

    data = read.table(each_file, header=F)
    TRA = data[data$V3=='TRA',]
    TRB = data[data$V3=='TRB',]

    p1 = ggplot(TRA,aes(x=reorder(TRA$V1,-TRA$V2),y=TRA$V2)) + 
      geom_bar(stat="identity",width = 0.5,fill =brewer.pal(7,"Set2")[1])+ # 修改柱条的宽度
      theme(axis.text.x=element_text(size=7,color="black",face="plain",vjust=0.5,hjust=0.5,angle=90)) +# Font face ("plain", "italic", "bold", "bold.italic") 
      geom_text(aes(label=TRA$V2,vjust=0.4,hjust=0),angle=90,size=2.5) +
      ylim(min(TRA$V2, 0)*1.1, max(TRA$V2)*1.2)+
      #scale_y_continuous(expand = c(0,0))+
      xlab("TRA genes") +ylab("Reads number") +
      ggtitle(paste(sample_name,"TRA",type_name,"reads number",sep=" "))

    p2 = ggplot(TRB,aes(x=reorder(TRB$V1,-TRB$V2),y=TRB$V2)) + 
      geom_bar(stat="identity",width = 0.5,fill=brewer.pal(7,"Set2")[2])+
      theme(axis.text.x=element_text(size=7,color="black",face="plain",vjust=0.5,hjust=0.5,angle=90)) +
      xlab("TRB genes") + ylab("Reads number")+
      ggtitle(paste(sample_name,"TRB",type_name,"reads number",sep=" "))+
      geom_text(aes(label=TRB$V2,vjust=0.4,hjust=0),angle=90,size=2.5)+
      ylim(min(TRB$V2, 0)*1.1, max(TRB$V2)*1.2)

    ggdraw() +
      draw_plot(p1, 0,.5, 1, .5) +
      draw_plot(p2, 0, 0, 1, .5) +
      draw_plot_label(c("A","B"),c(0,0),c(1,.5),size=15)

    output_name = paste(sample_name,type_name,"reads_number_barplot.jpg",sep=".")
    ggsave(paste(plot_dir,output_name,sep="/"))
}