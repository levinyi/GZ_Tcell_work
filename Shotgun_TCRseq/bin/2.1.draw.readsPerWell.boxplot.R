library(ggplot2)
library(reshape2)

args = commandArgs(T)

input_file = args[1]
print(paste("input:", input_file, sep=" "))
a = unlist(strsplit(input_file,"/"))
sample_name = unlist(strsplit(a[length(a)],"\\."))[2] # Total.G208E1.TRAB.wells.boole.matrix.Filtered.csv : extract G208E1
# print(sample_name)

data = read.table(input_file, sep = ",", header = T)
melted_data = melt(data) # transposition
# melted_data_Filtered = melted_data[which(melted_data$value>0),] # filtering 0
mdf = melted_data[melted_data$value>0,] # filtering

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
if(length(unique(mdf$variable))>24){
  axis_text_x = theme(axis.text.x = element_text(vjust = 0.5,hjust = 0.5, angle = 90))
}else{
  axis_text_x = theme(axis.text.x = element_text(vjust = 0.5,hjust = 0.5, angle = 0))
}

#####################################################
# draw plot 1
p = ggplot(mdf, aes(x=variable, y=value, fill=clone_type, color=clone_type)) + 
  geom_boxplot() + theme_bw()+
  xlab("Wells") +ylab("read count in wells")+
  theme(legend.title = element_blank()) +
  common_themes +
  axis_text_x
ggsave(paste(sample_name, "read_count.boxplot.png", sep = "_"), width = 10, height = 6)
print(paste("output: ", sample_name, "read_count.boxplot.png", sep="_"))
#################################################
# draw plot 2
min_data = aggregate(mdf[,4], by=list(variable=mdf$variable, clone_type=mdf$clone_type), FUN=min)
p1 = ggplot(min_data, aes(x=variable,y=x, group=clone_type, color=clone_type)) +
  geom_point() +
  geom_line() +
  theme_bw() +
  common_themes+
  axis_text_x +
  xlab("Wells") + ylab("The minimum read count in Wells")
ggsave(paste(sample_name, "Minimum_read_count_for_each_well.png", sep = "_"), width = 10, height = 6)
print(paste("output: ", sample_name, "Minimum_read_count_for_each_well.png", sep="_"))
