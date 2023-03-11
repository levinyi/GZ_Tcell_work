#' example: Rscript  /cygene/work/00.test/pipeline/10xGenomics/scripts/Infinitil_Plot_TCRs_on_UMAP.R integrated_seurat.rds tcr_id.txt
#' \input tcr_id.txt
#' \input integrated_seurat.rds
#' \output two pdfs.

library(Seurat)
library(ggplot2)
library(magrittr)
library(tibble)

# args = commandArgs(T)
# setwd("/cygene2/pipeline/10X/data/SA501001-TIL01/analysis/")
# setwd("/cygene2/pipeline/10X/data/SA501001-TIL02/analysis/")
setwd("/cygene2/pipeline/10X/data/CV501002/subsetTcells_re_analysis/analysis_w_BatchCor/")

# seurat_rds = args[1] # "integrated_seurat.rds"
# seurat_rds = "/cygene2/pipeline/10X/data/SA501001-TIL01/analysis/integrated_seurat.rds"
seurat_rds = "/cygene2/pipeline/10X/data/CV501002/subsetTcells_re_analysis/analysis_w_BatchCor/integrated_seurat.rds"

srt_obj = readRDS(seurat_rds)
# tcr_file = args[2]
# tcr_file = "/cygene2/pipeline/10X/data/SA501001-TIL01/analysis/SA511-TIL01_CD8_selected_TCR_PrintToUMAP.txt"
# tcr_file = "/cygene2/pipeline/10X/data/SA501001-TIL02/analysis/SA511-TIL02_tumor_reactive_TCRs.txt"
tcr_file = "/cygene2/pipeline/10X/data/CV501002/subsetTcells_re_analysis/analysis_w_BatchCor/selected_TCR_printToUmap.txt"
df = read.table(tcr_file, header=F)
tcr_id = df$V1
tcr_id
########################################
highlighted_barcode = rownames(srt_obj@meta.data[which(srt_obj@meta.data$clonotype_id %in% tcr_id),])
p1 = DimPlot(srt_obj, cells.highlight = highlighted_barcode)+ 
    scale_color_manual(
        values = c("#AEAAAA","#FF33CC"),
        breaks = c("Unselected","Group_1"),
        labels = c("FALSE","TRUE")
    )
p1
# ggsave("SA511-TIL01_tumor_reactive_TCRs.UMAPs.png", p1, device="png", width=6, height = 5)
ggsave("Selected_TCRs_PrintToUMAP.png", p1, device="png", width=6, height = 5)


#########################################
# for individual color
data = FetchData(srt_obj, vars=c("UMAP_1", "UMAP_2","clonotype_id"))

data$plot <- ifelse(!data$clonotype_id %in% tcr_id, "others", data$clonotype_id)
data$plot <- ifelse(grepl("group", data$plot), substring(data$plot,8), data$plot)

######### deal with clonotype name for TIL02 only
# 修改图例：
# special = c("clonotype1163","clonotype329")
data$plot <- ifelse(data$plot == "clonotype181", paste0(data$plot,"(against Hela,SiHA)"), data$plot)
data$plot <- ifelse(data$plot == "clonotype49", paste0(data$plot,"(against Hela, SiHA)"), data$plot)
data$plot <- ifelse(data$plot == "clonotype138", paste0(data$plot,"(against Hela)"), data$plot)

#########
library(ggrepel)
data_list = split(data, data$plot != 'others')
p2 = ggplot(data = data_list$`FALSE`,aes(x=UMAP_1,y=UMAP_2)) + geom_point(size=0.7,alpha=1, color='grey') +
    geom_point(data = data_list$`TRUE`, aes(x=UMAP_1,y=UMAP_2, color=plot),
               size=1.5) + 
    # geom_text_repel(data=data_list$`TRUE`,
    #                 label=data_list$`TRUE`$plot,
    #                 size = 2.5,
    #                 max.overlaps = length(rownames(data_list$`TRUE`)),
    #                 )+
    theme_classic() + 
    theme(
        legend.title = element_blank(),
        legend.position = "bottom"
    )+
    guides(
        color = guide_legend(
            ncol = 2
        )
    )
p2
# ggsave('SA511-TIL01_validated_CD4_TCR_plot_to_UMAPs.png',p2, device = "png",width = 6, height = 5)
ggsave('Selected_TCRs_PrintToUMAP.w.color.png',p2, device = "png",width = 6, height = 6)
# ggsave('SA511-TIL02_tumor_reactive_TCR_plot_w.text.png',p2, device = "png",width = 6, height = 5)
# ggsave('SA511-TIL02_tumor_reactive_TCR_plot_w.partial.text.png',p2, device = "png",width = 6, height = 5)
# ggsave('SA511-TIL02_tumor_reactive_TCR_plot_wo.text.png',p2, device = "png",width = 6, height = 5)

############################################################




##################################################
# validated_tcr = read.table("/cygene2/pipeline/10X/data/SA501001-TIL01/analysis/SA511-TIL01_validated_TCRs_PrintToUMAP.txt")
validated_tcr = read.table("/cygene2/pipeline/10X/data/SA501001-TIL02/analysis/SA511-TIL02_validated_TCRs_PrintToUMAP.txt")
head(validated_tcr)

data$validated <- ifelse(!data$clonotype_id %in% validated_tcr$V1, "others", data$clonotype_id)

for (i in data$validated){
    if (!i == "others"){
        data[which(data$validated == i),"validated"] <- validated_tcr[which(validated_tcr$V1==i),"V2"]
    }
}

data_list = split(data, data$validated != 'others')
p3 = ggplot(data = data_list$`FALSE`,aes(x=UMAP_1,y=UMAP_2)) + geom_point(size=0.7,alpha=1, color='#E7E6E6') +
    geom_point(data = data_list$`TRUE`, aes(x=UMAP_1,y=UMAP_2, color=validated),
               size=2) + 
    theme_classic() + 
    theme(
        legend.title = element_blank()
    ) + scale_color_manual(
        values = c("#C00000","#FFC000",'#5B9BD5','#AEAAAA'),
        breaks = c("Tumor-reactive","Self-reactive","Donor-reactive",'Negative'),
        labels = c("Tumor-reactive","Self-reactive","Donor-reactive",'Negative')
    )
p3
# ggsave('SA511-TIL01_validated_Tumor_reactive_TCRs.png',p3, device = "png",width = 6, height = 5)
ggsave('SA511-TIL02_validated_Tumor_reactive_TCRs.png',p3, device = "png",width = 6, height = 5)





#########################################################
metadata = c("CXCL13")
gene_expression_data = FetchData(srt_obj, vars=metadata)
head(gene_expression_data)
write.table(gene_expression_data %>% rownames_to_column('barcode'), file="SA511-TIL01_gene_umi_count.csv",sep=",",row.names = F,quote = T,col.names = T)
