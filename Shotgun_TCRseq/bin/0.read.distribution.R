library(ggplot2)
library(dplyr)

args = commandArgs(T)
# input_file1 = G208E1L2.barcode.total.rate.xls
# input_file2 = G208E1L1.barcode.total.rate.xls

# data combine
data1 = read.table(args[1])
data2 = read.table(args[2])

merged_data = rbind(data1, data2)

names(merged_data) = c("WellId", "barcode", "split_reads", "total_reads", "split_rate")

merged_data$well = substr(merged_data$WellId,5,7)  # 字符串截取
merged_data$clonotype = substr(merged_data$WellId, 3,3)

merged_data_2 = select(merged_data, well, clonotype, split_rate)


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
ggplot(merged_data_2, aes(x=well, y=split_rate, group=clonotype, color=clonotype)) + geom_point() + geom_line() +
	xlab("Well") + ylab("split rate") +
	ggtitle("split rate in each well") +
	common_themes

ggsave("split.rate.in.each.well.png", width=10, height=6)
print(paste("output:","split.rate.in.each.well.png",sep=" "))
