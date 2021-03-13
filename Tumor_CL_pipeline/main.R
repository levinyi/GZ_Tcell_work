library(tidyverse)
library(magrittr)
library(here)
here::i_am("main.R")
args = commandArgs(T)
source(here::here('src', 'Celligner_helpers.R'))
source(here::here('src', 'analysis_helpers.R'))
source(here::here('src', 'global_params.R'))


load_data <- function(data_dir, tumor_file = 'TCGA_mat.tsv', cell_line_file = 'CCLE_mat.csv', 
    annotation_file = 'Celligner_info.csv', hgnc_file = "hgnc_complete_set_7.24.2018.txt",
    RootPath_Tumor_file = "RootPath.Tumor.tpm.20210308.22.txt", RootPath_Tumor_info_file = "RootPath.Tumor.Info.20210308.22.csv",
    RootPath_CL_file = "RootPath.CellLine.tpm.20210308.2.txt", RootPath_CL_info_file = "RootPath.CellLine.Info.20210308.2.csv"){
    
    hgnc.complete.set <- data.table::fread(file.path(data_dir, hgnc_file)) %>% as.data.frame()
    
    TCGA_mat <-  readr::read_tsv(file.path(data_dir, tumor_file)) %>% as.data.frame() %>% tibble::column_to_rownames('Gene') %>% as.matrix() %>% t()
    CCLE_mat <-  readr::read_csv(file.path(data_dir, cell_line_file)) %>% as.data.frame() %>% tibble::column_to_rownames('X1') %>% as.matrix()
    colnames(CCLE_mat) <- stringr::str_match(colnames(CCLE_mat), '\\((.+)\\)')[,2]

    RTPT_tumor_mat = data.table::fread(file.path(data_dir, RootPath_Tumor_file)) %>% as.data.frame() %>% tibble::column_to_rownames('Geneid') %>% as.matrix() %>% t()
    RTPT_CL_mat = data.table::fread(file.path(data_dir, RootPath_CL_file)) %>% as.data.frame() %>% tibble::column_to_rownames('Geneid') %>% as.matrix() %>% t()
    RTPT_tumor_mat = log2(RTPT_tumor_mat+1) # log2(TPM+1)
    RTPT_CL_mat = log2(RTPT_CL_mat+1)

    common_genes <- intersect(colnames(TCGA_mat), hgnc.complete.set$symbol)
    common_genes <- intersect(common_genes, colnames(RTPT_tumor_mat))
    common_genes <- intersect(common_genes, colnames(RTPT_CL_mat))
    TCGA_mat <- TCGA_mat[,common_genes]
    RTPT_tumor_mat <- RTPT_tumor_mat[,common_genes]
    RTPT_CL_mat <- RTPT_CL_mat[,common_genes]

    hgnc.complete.set <- filter(hgnc.complete.set, symbol %in% common_genes)
    hgnc.complete.set <- hgnc.complete.set[!duplicated(hgnc.complete.set$symbol),]
    rownames(hgnc.complete.set) <- hgnc.complete.set$symbol
    hgnc.complete.set <- hgnc.complete.set[common_genes,]
    colnames(TCGA_mat) <- hgnc.complete.set$ensembl_gene_id
    colnames(RTPT_tumor_mat) <- hgnc.complete.set$ensembl_gene_id
    colnames(RTPT_CL_mat) <- hgnc.complete.set$ensembl_gene_id


    ann <- data.table::fread(file.path(data_dir, annotation_file)) %>% as.data.frame()
    RTPT_tumor_ann <- data.table::fread(file.path(data_dir, RootPath_Tumor_info_file)) %>% as.data.frame()
    RTPT_CL_ann <- data.table::fread(file.path(data_dir, RootPath_CL_info_file)) %>% as.data.frame()
    all_ann = rbind(ann,RTPT_tumor_ann,RTPT_CL_ann)

    TCGA_ann = dplyr::filter(all_ann, type=="tumor")
    CCLE_ann = dplyr::filter(all_ann, type=="CL")
    RTPT_tumor_ann = dplyr::filter(all_ann, type=="rootpath")
    RTPT_CL_ann = dplyr::filter(all_ann, type=="rootpath")

    func_genes <- dplyr::filter(hgnc.complete.set, !locus_group %in% c('non-coding RNA', 'pseudogene'))$ensembl_gene_id
    genes_used <- intersect(colnames(TCGA_mat), colnames(CCLE_mat))
    genes_used <- intersect(genes_used, colnames(RTPT_tumor_mat))
    genes_used <- intersect(genes_used, colnames(RTPT_CL_mat))
    genes_used <- intersect(genes_used, func_genes)

    TCGA_mat = TCGA_mat[,genes_used]
    CCLE_mat = CCLE_mat[,genes_used]
    RTPT_tumor_mat = RTPT_tumor_mat[,genes_used]
    RTPT_CL_mat = RTPT_CL_mat[,genes_used]


    return(list(TCGA_mat = TCGA_mat, TCGA_ann = TCGA_ann, CCLE_mat = CCLE_mat, CCLE_ann = CCLE_ann, 
      RTPT_tumor_mat = RTPT_tumor_mat, RTPT_tumor_ann = RTPT_tumor_ann, RTPT_CL_mat = RTPT_CL_mat, RTPT_CL_ann = RTPT_CL_ann))
}

create_Seurat_object <- function(exp_mat, ann, type = NULL) {
  seu_obj <- Seurat::CreateSeuratObject(t(exp_mat), min.cells = 0, min.features = 0, metadata = ann %>% magrittr::set_rownames(ann$sampleID))
  if(!is.null(type)) {
    seu_obj@meta.data$type <- type
  }
  # mean center the data, important for PCA
  seu_obj <- Seurat::ScaleData(seu_obj, features = rownames(Seurat::GetAssayData(seu_obj)), do.scale = F)
  seu_obj %<>% Seurat::RunPCA(assay='RNA', features = rownames(Seurat::GetAssayData(seu_obj)), npcs = global$n_PC_dims, verbose = F)
  seu_obj %<>% Seurat::RunUMAP(assay = 'RNA', dims = 1:global$n_PC_dims, reduction = 'pca', n.neighbors = global$umap_n_neighbors, min.dist =  global$umap_min_dist, metric = global$distance_metric, verbose=F)
  return(seu_obj)
}
cluster_data <- function(seu_obj) {
  seu_obj <- Seurat::FindNeighbors(seu_obj, reduction = 'pca', dims = 1:global$n_PC_dims, k.param = 20, force.recalc = TRUE, verbose = FALSE)
  seu_obj %<>% Seurat::FindClusters(reduction = 'pca', resolution = global$mod_clust_res)
  seu_obj@meta.data$cluster <- seu_obj@meta.data$seurat_clusters
  return(seu_obj)
}
#####################################
find_differentially_expressed_genes <- function(seu_obj) {
  n_clusts <- nlevels(seu_obj@meta.data$seurat_clusters)
  if (n_clusts > 2) {
    cur_DE_genes <- run_lm_stats_limma_group(t(Seurat::GetAssayData(seu_obj, assay='RNA', slot='scale.data')), seu_obj@meta.data %>% dplyr::select(seurat_clusters), limma_trend = TRUE) %>% dplyr::select(Gene, gene_stat = F_stat)
  } else if (n_clusts == 2) {
    cur_DE_genes <- run_lm_stats_limma(t(Seurat::GetAssayData(seu_obj, assay='RNA', slot='scale.data')), seu_obj@meta.data$cluster, limma_trend = TRUE) %>% dplyr::mutate(gene_stat = abs(t_stat)) %>% dplyr::select(Gene, gene_stat)
  } else {
    cur_DE_genes <- data.frame(Gene = colnames(seu_obj), gene_stat = NA)
  }
  return(cur_DE_genes)  
}
run_cPCA <- function(TCGA_obj, CCLE_obj, pc_dims = NULL) {
    cov_diff_eig <- run_cPCA_analysis(t(Seurat::GetAssayData(TCGA_obj, assay='RNA', slot='scale.data')), 
                                    t(Seurat::GetAssayData(CCLE_obj, assay='RNA', slot='scale.data')), 
                                    TCGA_obj@meta.data, CCLE_obj@meta.data, pc_dims=pc_dims)
 return(cov_diff_eig) 
}

# run mutual nearest neighbors batch correction
run_MNN <- function(CCLE_cor, TCGA_cor,  k1 = global$mnn_k_tumor, k2 = global$mnn_k_CL, ndist = global$mnn_ndist, subset_genes) {
  mnn_res <- modified_mnnCorrect(CCLE_cor, TCGA_cor, k1 = k1, k2 = k2, ndist = ndist, subset_genes = subset_genes)
  return(mnn_res)
}

# calculate the correlation between tumors and cell lines in the Celligner_aligned_data
calc_tumor_CL_cor <- function(Celligner_aligned_data, Celligner_info) {
  tumors_samples <- dplyr::filter(Celligner_info, type=='tumor')$sampleID
  cl_samples <- dplyr::filter(Celligner_info, type=='CL')$sampleID
  tumor_CL_cor <- cor(t(Celligner_aligned_data[tumors_samples,]), t(Celligner_aligned_data[cl_samples,]), use='pairwise')
  return(tumor_CL_cor)
}
calc_tumor_CL_cor_new <- function(Celligner_aligned_data, all_ann) {
    tumors_samples <- rbind(dplyr::filter(comb_ann, type=='tumor'), dplyr::filter(comb_ann, type=="rootpath"))$sampleID
    cl_samples <- dplyr::filter(comb_ann, type=='CL')$sampleID
    tumor_CL_cor <- cor(t(Celligner_aligned_data[tumors_samples,]), t(Celligner_aligned_data[cl_samples,]), use='pairwise')
    return(tumor_CL_cor)
}

calc_gene_stats <- function(TCGA_mat, CCLE_mat, hgnc_file="hgnc_complete_set_7.24.2018.txt") {
  common_genes <- intersect(colnames(TCGA_mat), colnames(CCLE_mat))

  hgnc.complete.set <- data.table::fread(file.path(data_dir, hgnc_file)) %>% as.data.frame()
  hgnc.complete.set <- hgnc.complete.set %>% 
    dplyr::select(Gene = ensembl_gene_id, Symbol = symbol) %>%
    filter(Gene %in% common_genes)
  hgnc.complete.set <- hgnc.complete.set[!duplicated(hgnc.complete.set$Gene),]
  rownames(hgnc.complete.set) <- hgnc.complete.set$Gene
  hgnc.complete.set <- hgnc.complete.set[common_genes,]  
  
  gene_stats <- data.frame(
    Tumor_SD = apply(TCGA_mat, 2, sd, na.rm=T),
    CCLE_SD = apply(CCLE_mat, 2, sd, na.rm=T),
    Tumor_mean = colMeans(TCGA_mat, na.rm=T),
    CCLE_mean = colMeans(CCLE_mat, na.rm=T),
    Gene = common_genes, stringsAsFactors = F) %>% dplyr::mutate(max_SD = pmax(Tumor_SD, CCLE_SD, na.rm=T)) #add avg and max SD per gene
  gene_stats <- left_join(hgnc.complete.set, gene_stats, by = "Gene")
  return(gene_stats)
}


##########################################
#############
#############       main
#############
##########################################

data_dir = "./"
RootPath_Tumor_file = args[1]
RootPath_info_file = args[2]
dat = load_data(data_dir, RootPath_Tumor_file=RootPath_Tumor_file, RootPath_info_file=RootPath_info_file)



comb_ann <- rbind(
  dat$TCGA_ann %>% dplyr::select(sampleID, lineage, subtype, `Primary/Metastasis`) %>% dplyr::mutate(type = 'tumor'),
  dat$CCLE_ann %>% dplyr::select(sampleID, lineage, subtype, `Primary/Metastasis`) %>% dplyr::mutate(type = 'CL')
)
# repeating paper:
# repeating Fig1:
source(here::here('src', 'Figure1.R'))
uncorrected_combined_type_plot = plot_uncorrected_data(dat$CCLE_mat, dat$TCGA_mat,comb_ann)
ggsave("Fig1.uncorrected_repeating.png")
# add rootpath data Fig1:
comb_ann <- rbind(
  dat$TCGA_ann %>% dplyr::select(sampleID, lineage, subtype, `Primary/Metastasis`) %>% dplyr::mutate(type = 'tumor'),
  dat$CCLE_ann %>% dplyr::select(sampleID, lineage, subtype, `Primary/Metastasis`) %>% dplyr::mutate(type = 'CL'),
  dat$RTPT_tumor_ann %>% dplyr::select(sampleID, lineage, subtype, `Primary/Metastasis`) %>% dplyr::mutate(type = 'rootpath'),
  dat$RTPT_CL_ann %>% dplyr::select(sampleID, lineage, subtype, `Primary/Metastasis`) %>% dplyr::mutate(type = 'rootpath')
)
uncorrected_combined_type_plot = plot_uncorrected_customized_data(dat$CCLE_mat, dat$TCGA_mat, comb_ann,dat$RTPT_tumor_mat, dat$RTPT_CL_mat)
ggsave("Fig1.uncorrected_combined_type_plot.with.my.tumor.png")
##################################################################################
# Fig2:
TCGA_obj <- create_Seurat_object(dat$TCGA_mat, dat$TCGA_ann, type='tumor')
CCLE_obj <- create_Seurat_object(dat$CCLE_mat, dat$CCLE_ann, type='CL')
TCGA_obj <- cluster_data(TCGA_obj)
CCLE_obj <- cluster_data(CCLE_obj)

TCGA_my_obj <- create_Seurat_object(rbind(dat$TCGA_mat, dat$RTPT_tumor_mat), rbind(dat$TCGA_ann, dat$RTPT_tumor_ann), type='tumor')
TCGA_my_obj <- cluster_data(TCGA_my_obj)
CCLE_my_obj <- create_Seurat_object(rbind(dat$CCLE_mat, dat$RTPT_CL_mat), rbind(dat$CCLE_ann, dat$RTPT_CL_ann), type='tumor') #### bug
CCLE_my_obj <- cluster_data(CCLE_my_obj)

tumor_DE_genes <- find_differentially_expressed_genes(TCGA_obj)
tumor_DE_genes_1 <- find_differentially_expressed_genes(TCGA_my_obj)
CL_DE_genes   <- find_differentially_expressed_genes(CCLE_obj)
CL_DE_genes_1 <- find_differentially_expressed_genes(CCLE_my_obj)

gene_stats   <- calc_gene_stats(dat$TCGA_mat, dat$CCLE_mat)
gene_stats_1 <- calc_gene_stats(rbind(dat$TCGA_mat,dat$RTPT_tumor_mat), rbind(dat$CCLE_mat,dat$RTPT_CL_mat))
##########
DE_genes <- full_join(tumor_DE_genes, CL_DE_genes, by = 'Gene', suffix = c('_tumor', '_CL')) %>%
    mutate(tumor_rank = dplyr::dense_rank(-gene_stat_tumor),
      CL_rank = dplyr::dense_rank(-gene_stat_CL),
      best_rank = pmin(tumor_rank, CL_rank, na.rm=T)) %>%
    dplyr::left_join(gene_stats, by = 'Gene')
# take genes that are ranked in the top 1000 from either dataset, used for finding mutual nearest neighbors
DE_gene_set <- DE_genes %>% dplyr::filter(best_rank < global$top_DE_genes_per) %>% .[['Gene']]
cov_diff_eig <- run_cPCA(TCGA_obj, CCLE_obj, global$fast_cPCA)
##########
DE_genes_1 <- full_join(tumor_DE_genes_1, CL_DE_genes_1, by = 'Gene', suffix = c('_tumor', '_CL')) %>%
    mutate(tumor_rank = dplyr::dense_rank(-gene_stat_tumor),
      CL_rank = dplyr::dense_rank(-gene_stat_CL),
      best_rank = pmin(tumor_rank, CL_rank, na.rm=T)) %>%
    dplyr::left_join(gene_stats_1, by = 'Gene')
DE_gene_set_1 <- DE_genes_1 %>% dplyr::filter(best_rank < global$top_DE_genes_per) %>% .[['Gene']]
cov_diff_eig_1 <- run_cPCA(TCGA_my_obj, CCLE_my_obj, global$fast_cPCA)
# how to compare differences between DE_genes and DE_genes_1 ??????


if(is.null(global$fast_cPCA)) {
    cur_vecs <- cov_diff_eig$vectors[, global$remove_cPCA_dims, drop = FALSE]
    cur_vecs_1 <- cov_diff_eig_1$vectors[, global$remove_cPCA_dims, drop = FALSE]
} else {
    cur_vecs <- cov_diff_eig$rotation[, global$remove_cPCA_dims, drop = FALSE]
    cur_vecs_1 <- cov_diff_eig_1$rotation[, global$remove_cPCA_dims, drop = FALSE]
}
# how to compare differences between cur_vecs and cur_vecs_1 ??????


TCGA_my_mat = rbind(dat$TCGA_mat, dat$RTPT_tumor_mat)
CCLE_my_mat = rbind(dat$CCLE_mat, dat$RTPT_CL_mat)
rownames(cur_vecs) <- colnames(dat$TCGA_mat)
rownames(cur_vecs_1) <- colnames(TCGA_my_mat)

TCGA_cor <- resid(lm(t(dat$TCGA_mat) ~ 0 + cur_vecs)) %>% t()
CCLE_cor <- resid(lm(t(dat$CCLE_mat) ~ 0 + cur_vecs)) %>% t()
TCGA_my_cor <- resid(lm(t(TCGA_my_mat) ~ 0 + cur_vecs_1)) %>% t()
CCLE_my_cor <- resid(lm(t(CCLE_my_mat) ~ 0 + cur_vecs_1)) %>% t()

mnn_res   <- run_MNN(CCLE_cor, TCGA_cor,        k1 = global$mnn_k_tumor, k2 = global$mnn_k_CL, ndist = global$mnn_ndist, subset_genes = DE_gene_set)
mnn_res_1 <- run_MNN(CCLE_my_cor, TCGA_my_cor,  k1 = global$mnn_k_tumor, k2 = global$mnn_k_CL, ndist = global$mnn_ndist, subset_genes = DE_gene_set_1)
# detected tied distances to neighbors, see ?'BiocNeighbors-ties'  ##### bug??????

combined_mat    <- rbind(mnn_res$corrected, CCLE_cor)
combined_my_mat <- rbind(mnn_res_1$corrected, CCLE_my_cor)
comb_obj <- create_Seurat_object(combined_mat, comb_ann)
comb_obj <- cluster_data(comb_obj)
# Warning: The following arguments are not used: reduction
comb_my_obj <- create_Seurat_object(combined_my_mat, comb_ann)
comb_my_obj <- cluster_data(comb_my_obj)

# fig2：
alignment =  Seurat::Embeddings(comb_obj, reduction = 'umap') %>% as.data.frame() %>% set_colnames(c('UMAP_1', 'UMAP_2')) %>% 
    rownames_to_column(var = 'sampleID') %>% left_join(comb_ann, by = 'sampleID')
alignment2 = Seurat::Embeddings(comb_my_obj, reduction = 'umap') %>% as.data.frame() %>% set_colnames(c('UMAP_1', 'UMAP_2')) %>% 
    rownames_to_column(var = 'sampleID') %>% left_join(comb_ann, by = 'sampleID')

source(here::here('src', 'Figure2.R'))
fig_plot = Celligner_alignment_plot(alignment)
ggsave("Fig2.raw.umap.png")

Seurat::DimPlot(comb_obj, reduction = 'umap',group.by = "orig.ident",colos=tissue_colors)
png("dimplot.png")

selected_data = alignment2[which(alignment2$type=='rootpath'),]
selected_data <- selected_data[!duplicated(selected_data$sampleID),]
fig_plot = Celligner_alignment_customized_plot(alignment2, selected_data)
ggsave("Fig2.all_comb.umap.png")
# save.image("main.Rdata")
##########################
# 来自于某blog
# 把UMAP的坐标放到metadata中：
# comb_obj <- AddMetaData(comb_obj, comb_obj@reductions$umap@cell.embeddings, col.name = colnames(comb_obj@reductions$umap@cell.embeddings))

# class_avg <- pbmc@meta.data %>%
#     group_by(RNA_snn_res.2) %>%
#     summarise(
#       UMAP_1 = median(UMAP_1),
#       UMAP_2 = median(UMAP_2)
#     )

# ggplot(pbmc@meta.data ,aes(x=UMAP_1,y=UMAP_2))+
#     geom_point(aes(color=RNA_snn_res.2))+
#     scale_color_manual(values = allcolour)+
#     geom_text(aes(label = RNA_snn_res.2), data = class_avg)+
#     theme(text=element_text(family="Arial",size=18)) +
#     theme(panel.background = element_rect(fill='white', colour='black'), 
#                      panel.grid=element_blank(), axis.title = element_text(color='black', family="Arial",size=18),
#                      axis.ticks.length = unit(0.4,"lines"), axis.ticks = element_line(color='black'), axis.ticks.margin = unit(0.6,"lines"),
#                      axis.line = element_line(colour = "black"), axis.title.x=element_text(colour='black', size=18), 
#                      axis.title.y=element_text(colour='black', size=18), axis.text=element_text(colour='black',size=18),
#                      legend.title=element_blank(), legend.text=element_text(family="Arial", size=18),
#                      legend.key=element_blank())+
#     theme(plot.title = element_text(size=22,colour = "black",face = "bold"))  + 
#     guides(colour = guide_legend(override.aes = list(size=5)))




# figure 3a

Celligner_aligned_data = Seurat::GetAssayData(comb_obj) %>% as.matrix() %>% t()
Celligner_aligned_my_data = Seurat::GetAssayData(comb_my_obj) %>% as.matrix() %>% t()
tumor_CL_cor = calc_tumor_CL_cor(Celligner_aligned_data, comb_ann)
tumor_CL_my_cor = calc_tumor_CL_cor(Celligner_aligned_my_data, comb_ann)
write.table(t(tumor_CL_cor) %>% as.data.frame() %>% rownames_to_column("ccle_name"), 
    file="tumor_CL_cor.transit.csv", sep=",", row.names = FALSE, quote = FALSE)
write.table(t(tumor_CL_my_cor) %>% as.data.frame() %>% rownames_to_column("ccle_name"), 
    file="tumor_CL_my_cor.transit.csv", sep=",", row.names = FALSE, quote = FALSE)

# fig3
get_cell_line_tumor_class <- function(tumor_CL_cor, alignment) {
  cl_tumor_classes <- apply(tumor_CL_cor, 2, function(x) cell_line_tumor_class(x, tumor_CL_cor, alignment)) %>% as.character()
  names(cl_tumor_classes) <- colnames(tumor_CL_cor)
  return(cl_tumor_classes)
}
cl_tumor_classes = get_cell_line_tumor_class(tumor_CL_cor, alignment)
cell_line_tumor_class_plot(cl_tumor_classes, alignment, tumor_CL_cor, "Fig3a.myheatmap.png")
cl_tumor_classes2 = get_cell_line_tumor_class(tumor_CL_my_cor, alignment2)
cell_line_tumor_class_plot(cl_tumor_classes2, alignment2, tumor_CL_my_cor, "Fig3a.myheatmap.my.png")
# Fig 3a end

# Fig 3b start
fig3b_plot = cell_line_tumor_distance_distribution(alignment, tumor_CL_cor)
ggsave("Fig3b.tumor_dist_spread.png")

fig3b_plot = cell_line_tumor_distance_distribution(alignment2, tumor_CL_my_cor)
ggsave("Fig3b.tumor_dist_spread.my.png")
