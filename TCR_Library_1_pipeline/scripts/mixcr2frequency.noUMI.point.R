library(ggplot2)
library(RColorBrewer)
args = commandArgs(T)

# input file: GXXX.filtered2reads.freq.xls
a = unlist(strsplit(args[1],"/"))
sample_name = unlist(strsplit(a[length(a)], "\\."))[1]

theme_x_axis = theme(
  axis.text.x = element_text(size = 15),
  axis.text.y = element_text(size = 15),
)
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

data = read.table(args[1], header = T)
max_length = length(data$Frequency_Log10)
data$id = 1:max_length

p = ggplot(data, aes(x = id, y = Frequency_Log10)) + geom_point(color = brewer.pal(7, "Set1")[2]) +
	common_themes + theme_x_axis + 
	xlab("clonotypes") +
	ylab("Frequency log10") +
	ggtitle(sample_name)
ggsave(paste(sample_name, "filtered2reads.frequency.png", sep="."), width = 10, height = 6)
