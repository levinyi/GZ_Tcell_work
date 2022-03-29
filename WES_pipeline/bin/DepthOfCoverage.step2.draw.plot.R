library(tidyr)
library(dplyr)
library(ggplot2)
library(hrbrthemes) # for theme_ipsum()

# args = c("HC22_TU_CL1-P0-gDNA.Tumor.coverage.depth.rate.xls")
# chr     chr_len total_depth     coverage_leng   rate1   rate2
# chr1    4046489 234630948       3795965 57.983834380867 0.938088550345744
# chr10   1488150 86518107        1464131 58.138028424554 0.983859825958405

args = commandArgs(T)


theme_x_axis = theme(
  axis.text.x = element_text(size = 20),
  axis.text.y = element_text(size = 20),
)

common_themes = theme(plot.title = element_text(size = 20, hjust = 0.5),
                      panel.grid.major = element_line(colour = NA), # 去掉网格线
                      panel.background = element_blank(), # 去掉背景
                      panel.grid.major.x = element_line(),  # 横向网格线linetype=2,xuxian
                      panel.grid.major.y = element_line(), # 纵向网格线
                      panel.border = element_rect(fill=NA), # 边框
                      axis.line = element_line(), # 坐标轴线
                      axis.title = element_text(size = 20),
                      legend.title = element_blank(),
                      legend.position = c(0.95,0.95), # bottom, right, 'left'
                      legend.background = element_rect(colour = NA, fill = NA),
)
for (each_file in args){
  output = dirname(each_file)
  file_name = basename(each_file)
  title_name = paste(unlist(strsplit(file_name, "\\."))[1:2], collapse=".")
  print(title_name)
  data = read.table(each_file, header=T)
  # print(data)
  max_data = max(data$rate1)
  print(paste("max_data=",max_data,sep = ""))
  if (max_data>=100){
    t = floor(max_data*1.1)
  }else{
    t = 100
  }
  print(paste("t=",t,sep=""))
  p = ggplot(data) + 
      geom_bar(aes(x=chr, y=rate1), stat = 'identity',fill = "#8CC53E", width = 0.7)+
      geom_line( aes(x=chr, y=rate2*t/1),group=1) +
  		geom_point(aes(x=chr, y=rate2*t/1))+
  		xlab("")+
  		theme_bw() +theme_x_axis + common_themes +
  		theme(axis.text.x = element_text(angle = 90, hjust=0.1, vjust=0.5))+
      scale_y_continuous(name=expression("Mean Depth"),
                         limits = c(0,t),breaks = seq(0,t,20),
                         sec.axis = sec_axis(~./t, name = "Proportion of Covered Bases"))
  ggsave(filename = paste(paste(output,title_name,sep="/"), "DepthOfCoverage.png", sep="."), plot = p, width=9, height=7, path = "./" )
  print(paste(title_name,"DepthOfCoverage.png", sep="."))
}
print("all done")


# ggplot(data,aes(x=chr)) +
#   geom_bar( aes(y=rate1), stat = 'identity',size=.1, fill = "#8CC53E", width = 0.7)+
#   geom_line( aes(y=rate2*t), group=1, size=2, color="#69b3a2") +
#   geom_point(aes(x=chr, y=rate2*t)) +
#   xlab("")+
#   theme_ipsum_pub()+
#   theme_bw() +theme_x_axis + common_themes +
#   theme(axis.text.x = element_text(angle = 90, hjust=0.1, vjust=0.5))+
#   scale_y_continuous(
#     # Features of the first axis
#     name=expression("Mean Depth"),
#     limits = c(0,t),breaks = seq(0,t,20),
#     
#     # Add a second axis and specify its features
#     sec.axis = sec_axis(~./t, name = "Proportion of Covered Bases")) + 
#   theme(
#     axis.title.y = element_text(color = "#8CC53E", size=13),
#     axis.title.y.right = element_text(color = "#69b3a2", size=13)
#   ) +
#   ggtitle("DepthOfCoverage")
# library('hrbrthemes')
# install.packages("hrbrthemes")
