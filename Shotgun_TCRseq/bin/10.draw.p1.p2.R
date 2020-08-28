library(ggplot2)
library(RColorBrewer)

args = commandArgs(T)

theme_x_axis = theme(
  axis.text.x = element_text(size = 15),
  axis.text.y = element_text(size = 15),
)
common_themes = theme(
	plot.title = element_text(size = 17, hjust = 0.5),
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

sample_name = args[1]
rawdata = read.table(args[2], header = T) # "Total.pairs.FromA2B.threshold.50.add.xls"

data_p1 = rawdata[order(rawdata$p1),]
data_p2 = rawdata[order(rawdata$p2),]

data_p1$id = 1:length(data_p1$TRA_clone)
data_p2$id = 1:length(data_p2$TRA_clone)


ggplot(data_p1, aes(x=id, y=p1)) + geom_point(color=brewer.pal(7,"Set1")[2]) +
  common_themes +
  xlab("clonotype") + ylab("p1 value") +
  scale_x_continuous(limits = c(0,20500)) +
  scale_y_continuous(limits = c(0,0.2)) +
  theme_x_axis +ggtitle(sample_name)
ggsave(paste(sample_name, "p1.png", sep="."), width=10, height=6)

ggplot(data_p2, aes(x=id, y=p2)) + geom_point(color=brewer.pal(7,"Set1")[2]) +
  common_themes +
  xlab("clonotype") + ylab("p2 value") +
  scale_x_continuous(limits = c(0,25000)) +
  scale_y_continuous(limits = c(0,1.1), breaks = seq(0,1.25,0.1)) +
  theme_x_axis + ggtitle(sample_name)
ggsave(paste(sample_name, "p2.png", sep="."), width=10, height=6)
