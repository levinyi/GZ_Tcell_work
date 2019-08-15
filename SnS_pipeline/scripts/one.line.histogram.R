library(ggplot2)

args = commandArgs(T)

data = read.table(args,header=F)

ggplot(data,aes(x=data$V1)) + geom_histogram(binwidth=4,stat="count")
ggsave("histogram.jpg")
