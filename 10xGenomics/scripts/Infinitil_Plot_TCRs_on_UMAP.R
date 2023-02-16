#' example: Rscript  /cygene/work/00.test/pipeline/10xGenomics/scripts/Infinitil_Plot_TCRs_on_UMAP.R integrated_seurat.rds tcr_id.txt
#' \input tcr_id.txt
#' \input integrated_seurat.rds
#' \output two pdfs.

library(Seurat)
library(ggplot2)
library(magrittr)

args = commandArgs(T)

seurat_rds = args[1] # "integrated_seurat.rds"
srt_obj = readRDS(seurat_rds)

tcr_file = args[2]
tcr_id = read.table(tcr_file, header=F)
tcr_id = tcr_id$V1
# extract barcode
highlighted_barcode <- rownames(srt_obj@meta.data[srt_obj@meta.data$clonotype_id %in% tcr_id,])
print(head(highlighted_barcode))
## 
p <- DimPlot(srt_obj, cells.highlight = highlighted_barcode)
p <- p + scale_color_manual(
		values=c("grey50","red"),
		breaks=c("Unselected","Group_1"),
		labels=c("Others","Selected")
	)
ggsave("Plot.Selected.TCRs.on.UMAP.pdf",plot=p, device="pdf",width=6.2, height=5)
ggsave("Plot.Selected.TCRs.on.UMAP.png",plot=p, device="png",width=6.2, height=5)

#########################################
# for individual color
data = FetchData(srt_obj,vars=c("UMAP_1", "UMAP_2", "clonotype_id"))
# umap_coords <- as.data.frame(srt_obj@reductions$umap@cell.embeddings)

# add column to store tcr info.
data$tcr_clonotype <- "others"
for (each in tcr_id) {
	data$tcr_clonotype[which(data$clonotype_id == each)] <- each
}

print(head(data))
# split selected data and Unselected data. give ggplot2 two layers.
selected_data = data %>% dplyr::filter(tcr_clonotype != "others")
unselect_data = data %>% dplyr::filter(tcr_clonotype == "others")

# ggplot2
p2 <- ggplot(unselect_data, aes(x=UMAP_1,y=UMAP_2)) +
	geom_point(size=0.7, alpha=1, color="grey") +
	geom_point(data=selected_data, aes(x=UMAP_1,y=UMAP_2, color=tcr_clonotype),size=1) +
	theme_classic() + theme(legend.title=element_blank())

# convert group1_clonotype80 to clonotype80: 
short_name = c()
for (each in tcr_id) {
	if (!startsWith(each, "clonotype")){
		print(paste("this clonotype is not starts with clonotype", each, sep=":"))
		name = substring(each, 8)
		short_name = c(short_name, name)
	}
}
print(short_name)
p2 <- p2 + scale_color_discrete(breaks=tcr_id, labels=short_name)
ggsave("Plot.Selected.TCRs.on.UMAP.sep.cor.pdf", plot=p2, device="pdf", width=7.5,height=5)
ggsave("Plot.Selected.TCRs.on.UMAP.sep.cor.png", plot=p2, device="png", width=7.5,height=5)
print("done")
# library(ggrepel)
# p+ geom_text_repel(data = selected_data)
