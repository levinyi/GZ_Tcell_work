library(ggplot2)

args = commandArgs()

# input_file1 = G208E1L2.barcode.total.rate.xls
# input_file2 = G208E1L1.barcode.total.rate.xls

# data combine
data1 = read.table(args[1], sep="\t", header = F)
data2 = read.table(args[2], sep="\t", header = F)

data1$name = 
