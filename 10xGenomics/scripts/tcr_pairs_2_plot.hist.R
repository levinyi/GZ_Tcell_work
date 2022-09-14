library(ggplot2)
args = commandArgs(T)

# data from : tcr_table_2_pairs_format.py

df <- read.table(args[1])
sample_name <- args[2]
print(df)

p = ggplot(df, aes(x=V1,y=V2)) + geom_bar(stat="identity", fill="darkblue",) +
  geom_text(aes(label=V2, vjust=-1 ))+ ggtitle(sample_name) +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5)) + xlab("TCR types") +
  ylab("Frequency")
ggsave(paste(sample_name, "TCR_pairs.png", sep="."), p)
