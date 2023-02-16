library(Seurat)
args = commandArgs(T)
seurat_rds = args[1] # "integrated_seurat.rds"

srt_obj = readRDS(seurat_rds)

tcr_file = args[2]
tcr_id = read.table(tcr_file,header=F)


# extract barcode
highlighted_barcode <- rownames(srt_obj@meta.data[srt_obj@meta.data$clonotype_id %in% tcr_id,])

## 
p <- DimPlot(srt_obj, cells.highlight = highlighted_barcode)
p <- p + scale_color_manual(
		values=c("grey50","red"),
		breaks=c("Unselected","Group_1"),
		labels=c("Others","Selected")
	)
ggsave("test_Plot.Selected.TCRs.on.UMAP.pdf",plot=p, device="pdf",width=5, height=5)
#########################################
# for individual color
data = FetchData(srt_obj,vars=("UMAP_1","UMAP_2","clonotype_id"))
# umap_coords <- as.data.frame(srt_obj@reductions$umap@cell.embeddings)

# add column to store tcr info.
data$tcr_clonotype <- "others"
for (each in tcr_id) {
	data$tcr_clonotype[which(data$clonotype_id == each)] <- each
}

# split selected data and Unselected data. give ggplot2 two layers.
selected_data = data %>% filter(tcr_clonotype != "others")
unselect_data = data %>% filter(tcr_clonotype == "others")

# ggplot2
p2 <- ggplot(unselect_data, aes(x=UMAP_1,y=UMAP_2)) +
	geom_point(size=0.7, alpha=1, color="grey") +
	geom_point(data=selected_data, aes(x=UMAP_1,y=UMAP_2, color=tcr_clonotype),size=1) +
	theme_classic() + theme(legend.title=element_blank())
ggsave("test_Plot.Selected.TCRs.on.UMAP.sep.cor.pdf", plot=p2, device="pdf", width=8,height=5)

