library(tidyr)
library(ggplot2)
args = commandArgs(T)


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

for (each in args){
	path_list = unlist(strsplit(each_file,"/"))
	file_name = path_list[legth(path_list)]
	title_name = unlist(strsplit(file_ame, "\\."))[1:3]
  data = read.csv(each_file, header=T)
  data3 = data %>% separate(Locus,c("chr","pos"),":") %>%
  # separate(chr,c("chr1","chr"),"chr") %>%  
  group_by(chr) %>%
  summarise(chr_len = n(), total_depth = sum(Total_Depth),
            coverage_leng = sum(Total_Depth!=0)) %>%
  mutate(rate1 = total_depth/chr_len, rate2=coverage_leng/chr_len)

p = ggplot(data3) +
  geom_line( aes(x=chr, y=rate2*500/1),group=1) +
  geom_point(aes(x=chr, y=rate2*500/1)) +
  geom_bar(aes(x=chr, y=rate1), stat = 'identity',fill = "#8CC53E", width = 0.7) +
  xlab("")+
  scale_y_continuous(name=expression("Mean Depth"), limits = c(0, 500),
                     sec.axis = sec_axis(~./500, name = "Proportion of Covered Bases")) +
  theme_bw() +theme_x_axis + common_themes

ggsave(filename = paste(title_name,"coverage.png",sep="."),plot = p,height = 8,width = 10,)
write.table(data3, file=paste(title_name,"coverage.depth.rate.xls",sep="."), sep="\t", row.names = FALSE, quote=FALSE)
}
