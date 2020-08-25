library(ggplot2)
library(reshape2)
args = commandArgs(T)

# input file: Raw.reads.in.wells.mixcr.xls

data = read.table(args[1], sep = ",", header = T)
melt_data = melt(data)

TRA_data = melt_data[melt_data$cloneType=="TRA",]
TRA_data = TRA_data[TRA_data$value>=18, ]

TRB_data = melt_data[melt_data$cloneType=="TRB",]
TRB_data = TRB_data[TRB_data$value>=18,]
# head(TRB_data)

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

p1 = ggplot(TRA_data, aes(x=reorder(TRA_data$cloneId,-TRA_data$value), y=log(TRA_data$value), color=TRA_data$variable)) + geom_point() +geom_line()+
  theme(axis.text.x = element_text(vjust = 0.5,hjust = 0.5, angle = 90)) +
  ggtitle("Reads number in Each Well/clonotype (TRA)") +
  xlab("clone id") + ylab("reads number")+
  common_themes
ggsave("Reads.number.in.each.Well.clonotype.TRA.png", width=10, height=6)

p2 = ggplot(TRB_data, aes(x=reorder(TRB_data$cloneId,-TRB_data$value), y=log(TRB_data$value), color=TRB_data$variable)) + geom_point() +geom_line()+
  theme(axis.text.x = element_text(vjust = 0.5,hjust = 0.5, angle = 90)) +
  ggtitle("Reads number in Each Well/clonotype (TRB)") +
  xlab("clone id") + ylab("reads number")+
  common_themes
ggsave("Reads.number.in.each.Well.clonotype.TRB.png", width=10, height=6)
