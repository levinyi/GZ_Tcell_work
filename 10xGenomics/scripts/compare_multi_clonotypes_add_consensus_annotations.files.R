library(tidyverse)
library(VennDiagram)
library(RColorBrewer)
file1 = "C:\\Users\\cy006\\Nutstore\\.nutstore_c2hpeWlAcm9vdHBhdGhneC5jb20=\\DuShiYi\\P0000-blackbird\\2103-CR001\\CR001004\\CR001004_CD4_CD8_scRNA_TCRseq\\CR001004_G179E3L1_TCRseq_IMGT\\clonotypes_add_consensus_annotations.csv"
file2 = "C:\\Users\\cy006\\Nutstore\\.nutstore_c2hpeWlAcm9vdHBhdGhneC5jb20=\\DuShiYi\\P0000-blackbird\\2103-CR001\\CR001004\\CR001004_G179E1L1_TCRseq_IMGT\\clonotypes_add_consensus_annotations.csv"

data1 = read.csv(file1, sep=",", header=T)
data2 = read.csv(file2, sep=",", header=T)

pair_identify <- function(data){
  data = data %>% 
    tidyr::unite(clonotype_tra_id, c(v_gene, cdr3, j_gene, ), sep = "", remove = FALSE) %>% 
    tidyr::unite(clonotype_trb_id, c(v_gene.1, cdr3, j_gene.1,), sep = "", remove = FALSE)
  data = data %>% unite(clonotype_pair_id, c(clonotype_tra_id, clonotype_trb_id), sep = "_")
  return(data)
}
data1 = pair_identify(data1)
data2 = pair_identify(data2)

# draw venn diagram
file_name1 = "G179E3L1-PBMC-CD8"
file_name2 = "G179E1L1-TILs-CD3T"
venn.diagram(list(A=data1$clonotype_pair_id, A=data2$clonotype_pair_id),
             resolution = 300, imagetype = "png",
             col = RColorBrewer::brewer.pal(7, "Dark2")[1:2],
             cex = 4, # size
             label.col=c("black"),
             lwd = 4,
             lty = 1,
             # for category:
             cat.fontface="bold", # font type.
             category.names = c(file_name1, file_name2),
             cat.cex = 3.5,
             # fill = c("blue","red"), # no fill is more beautiful.
             cat.col = RColorBrewer::brewer.pal(7, "Dark2")[1:2],
             area.vector = TRUE,
             # for main title:
             main="", main.cex = 2,main.col = "black", main.fontface=2,
             sub.fontfamily = "serif", main.just = 0,
             # for subtitle:
             sub = "", sub.cex = 4, sub.col = "black", sub.fontface = 2,
             filename ="VennDiagram.pairwise.png")
##################

data_merged = dplyr::full_join(data1, data2, by = "clonotype_pair_id") %>% 
  select(clonotype_pair_id,clonotype_id.x,clonotype_id.y,frequency.x,frequency.y) %>% 
  drop_na() %>% dplyr::rename(
    paste("clonotype_id", file_name1) = clonotype_id.x,
    paste("clonotype_id", file_name2) = clonotype_id.y
  )
setwd("C:\\Users\\cy006\\Nutstore\\.nutstore_c2hpeWlAcm9vdHBhdGhneC5jb20=\\DuShiYi\\P0000-blackbird\\2103-CR001\\CR001004")
write.table(data_merged, file = "hhh.csv", sep = ",", row.names = FALSE)

