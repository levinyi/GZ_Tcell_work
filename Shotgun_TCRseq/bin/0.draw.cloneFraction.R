library(ggplot2)
library(RColorBrewer)

args = commandArgs(T)
# input file : G208E1L2.mixcr.out.clonotypes.TRA.txt
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
theme_x_axis = theme(
  axis.text.x = element_text(size = 15),
  axis.text.y = element_text(size = 15),
)

data = read.table(args[1], header = T, sep = "\t")

a = unlist(strsplit(args[1], "/"))
sample_name = unlist(strsplit(a[length(a)],"\\."))[1]  # G208E1L1
clonetype = unlist(strsplit(a[length(a)], "\\."))[5]  # TRA

max_length = length(data$cloneId)
data$id = 1:max_length

p = ggplot(data)+ geom_point(aes(x=id,y=cloneFraction,), color=brewer.pal(7,"Set1")[2]) +
  common_themes + 
  scale_y_log10(limits = c(1e-7, 1),) +
  scale_x_continuous(limits = c(0,max_length+200), breaks = seq(0,max_length+200,200))+
  xlab(paste(clonetype, "ranking", sep=" ")) +
  ggtitle(paste(sample_name, clonetype, "cloneFraction", sep=" ")) + theme_x_axis

output_name = paste(sample_name, clonetype, "cloneFraction.png", sep=".")

ggsave(output_name)
print(paste("output:", output_name, sep=" "))
