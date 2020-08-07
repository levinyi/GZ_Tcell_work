library(ggplot2)
args = commandArgs(T)

for (each_file in args){
    print(each_file)
    data = read.table(each_file, sep = ",", header = T)
    tra_data = subset(data, clone_type == 'TRA')[,c(-1,-2)]
    trb_data = subset(data, clone_type == 'TRB')[,c(-1,-2)]

    tra_data_sum = apply(tra_data, 2, sum)
    trb_data_sum = apply(trb_data, 2, sum)

    tra_data_sum_t = as.data.frame(tra_data_sum)
    tra_data_sum_t$barcode = row.names(tra_data_sum_t)
    tra_data_sum_t$group = c("TRA")

    trb_data_sum_t = as.data.frame(trb_data_sum)
    trb_data_sum_t$barcode = row.names(trb_data_sum_t)
    trb_data_sum_t$group = c("TRB")
    names(tra_data_sum_t) = c("data_sum","barcode","groups")
    names(trb_data_sum_t) = c("data_sum","barcode","groups")

    total_data = rbind(tra_data_sum_t, trb_data_sum_t)
    common_themes = theme(plot.title = element_text(hjust = 0.5),
                panel.grid.major = element_line(colour = NA), # 去掉网格线
                panel.background = element_blank(), # 去掉背景
                panel.grid.major.x = element_line(),  # 横向网格线linetype=2,xuxian
                panel.grid.major.y = element_line(), # 纵向网格线
                panel.border = element_rect(fill=NA), # 边框
                axis.line = element_line(), # 坐标轴线
                legend.title = element_blank()
          )

    if (endsWith(each_file, ".wells.count.matrix.csv")){
        p = ggplot(total_data, aes(x=barcode, y=data_sum, group=groups,color=groups)) +geom_point() + geom_line() + theme_bw() +
          xlab("Well Position") + ylab("Total Raw Read Count") +
          ggtitle("Total Raw read count for Each Well") +
          common_themes +
          scale_fill_discrete(breaks=c("TRA","TRB"), labels=c("TRA-Raw", "TRB-Raw"))
        ggsave( "read_count_for_each_well.Raw.jpg", width=10, height=10)
        print("output : read_count_for_each_well.Raw.jpg")
    }else if (endsWith(each_file, ".wells.boole.matrix.csv")){
        p = ggplot(total_data, aes(x=barcode,y=data_sum, group=groups, color=groups)) + geom_point() + geom_line() + theme_bw() +
            xlab("Well Position") + ylab("TCR Raw Count") +
            ggtitle("Number of TCRs detected for each well") +
            common_themes +
            scale_y_continuous(breaks = seq(0,100,5)) +
            scale_fill_discrete(breaks=c("TRA","TRB"), labels=c("TRA-Raw", "TRB-Raw"))
        ggsave( "TCRs_count_for_each_well.Raw.jpg", width=10, height=10)
        print("TCRs_count_for_each_well.Raw.jpg")
    }else if (endsWith(each_file, ".wells.count.matrix.Filtered.csv")){
        p = ggplot(total_data, aes(x=barcode, y=data_sum, group=groups,color=groups)) + geom_point() + geom_line() + theme_bw() +
          xlab("Well Position") + ylab("Total Filtered Read Count") +
          ggtitle("Total Filtered read count for Each Well") +
          common_themes +
          scale_fill_discrete(breaks=c("TRA","TRB"), labels=c("TRA-Filtered", "TRB-Filtered"))
        ggsave( "read_count_for_each_well.Filtered.jpg", width=10, height=10)
        print("read_count_for_each_well.Filtered.jpg")
    }else if (endsWith(each_file, ".wells.boole.matrix.Filtered.csv")){
        p = ggplot(total_data, aes(x=barcode, y=data_sum, group=groups, color=groups)) + geom_point() + geom_line() + theme_bw() +
          xlab("Well Position") + ylab("TCR Filtered Count") +
          ggtitle("Number of TCRs detected for each well") +
          common_themes +
          scale_y_continuous(breaks = seq(0,100,5)) +
          scale_fill_discrete(breaks=c("TRA","TRB"), labels=c("TRA-Filtered", "TRB-Filtered"))
        ggsave( "TCRs_count_for_each_well.Filtered.jpg", width=10, height=10)
        print("TCRs_count_for_each_well.Filtered.jpg")
    }else{
        exit("your input file name are not correct.\nError: Incorrect File Name.\n")
    }
}

