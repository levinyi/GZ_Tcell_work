setwd("/cygene/work/00.test/pipeline/Tumor_CL_pipeline")
library(here)
here()
library(magrittr)
library(tidyverse)
source(here::here('src', 'Celligner_helpers.R'))
source(here::here('src', 'analysis_helpers.R'))
source(here::here('src', 'global_params.R'))



load_data <- function(data_dir, tumor_file = 'database/TCGA_mat.tsv', cell_line_file = 'database/CCLE_mat.csv', 
                      annotation_file = 'database/Celligner_info.csv', hgnc_file = "database/hgnc_complete_set_7.24.2018.txt") {
  hgnc.complete.set <- data.table::fread(file.path(data_dir, hgnc_file)) %>% as.data.frame()
  
  
  TCGA_mat <-  readr::read_tsv(file.path(data_dir, tumor_file)) %>% as.data.frame() %>% tibble::column_to_rownames('Gene') %>% as.matrix() %>% t()
  
  common_genes <- intersect(colnames(TCGA_mat), hgnc.complete.set$symbol)
  TCGA_mat <- TCGA_mat[,common_genes]
  hgnc.complete.set <- filter(hgnc.complete.set, symbol %in% common_genes)
  hgnc.complete.set <- hgnc.complete.set[!duplicated(hgnc.complete.set$symbol),]
  rownames(hgnc.complete.set) <- hgnc.complete.set$symbol
  hgnc.complete.set <- hgnc.complete.set[common_genes,]
  colnames(TCGA_mat) <- hgnc.complete.set$ensembl_gene_id
  
  CCLE_mat <-  readr::read_csv(file.path(data_dir, cell_line_file)) %>% as.data.frame() %>% tibble::column_to_rownames('X1') %>% as.matrix()
  
  colnames(CCLE_mat) <- stringr::str_match(colnames(CCLE_mat), '\\((.+)\\)')[,2]
  
  if(is.null(annotation_file) | !file.exists(file.path(data_dir, annotation_file))) {
    ann <- data.frame(sampleID = c(rownames(TCGA_mat), rownames(CCLE_mat)),
                      lineage = NA,
                      subtype = NA,
                      type = c(rep('tumor', nrow(TCGA_mat)), rep('CL', nrow(CCLE_mat))))
    ann$`Primary/Metastasis` <- NA
  } else {
    ann <- data.table::fread(file.path(data_dir, annotation_file)) %>% as.data.frame()
    if('UMAP_1' %in% colnames(ann)) {
      ann <- ann %>% 
        dplyr::select(-UMAP_1)
    }
    if('UMAP_2' %in% colnames(ann)) {
      ann <- ann %>% 
        dplyr::select(-UMAP_2)
    }
    if('cluster' %in% colnames(ann)) {
      ann <- ann %>% 
        dplyr::select(-cluster)
    }
  }
  
  TCGA_ann <- dplyr::filter(ann, type=='tumor')
  CCLE_ann <- dplyr::filter(ann, type=='CL')
  
  func_genes <- dplyr::filter(hgnc.complete.set, !locus_group %in% c('non-coding RNA', 'pseudogene'))$ensembl_gene_id
  genes_used <- intersect(colnames(TCGA_mat), colnames(CCLE_mat))
  genes_used <- intersect(genes_used, func_genes)
  
  TCGA_mat <- TCGA_mat[,genes_used]
  CCLE_mat <- CCLE_mat[,genes_used]
  
  return(list(TCGA_mat = TCGA_mat, TCGA_ann = TCGA_ann, CCLE_mat = CCLE_mat, CCLE_ann = CCLE_ann))
}

load_my_data <- function(data_dir, hgnc_file = "database/hgnc_complete_set_7.24.2018.txt",
    RootPath_Tumor_file = "database/RootPath.Tumor.tpm.20210801.23.txt", RootPath_Tumor_info_file = "database/RootPath.Tumor.Info.20210801.23.csv",
    RootPath_CL_file = "database/RootPath.CellLine.tpm.20210308.2.txt", RootPath_CL_info_file = "database/RootPath.CellLine.Info.20210308.2.csv"){
    
    hgnc.complete.set <- data.table::fread(file.path(data_dir, hgnc_file)) %>% as.data.frame()
    RTPT_tumor_mat = data.table::fread(file.path(data_dir, RootPath_Tumor_file)) %>% as.data.frame() %>% tibble::column_to_rownames('Geneid') %>% as.matrix() %>% t()
    RTPT_CL_mat = data.table::fread(file.path(data_dir, RootPath_CL_file)) %>% as.data.frame() %>% tibble::column_to_rownames('Geneid') %>% as.matrix() %>% t()
    
    RTPT_tumor_mat = log2(RTPT_tumor_mat+1) # log2(TPM+1)
    RTPT_CL_mat = log2(RTPT_CL_mat+1)
    
    # 将RTPT_tumor_mat的列名的gene名称转化为ensembl_gene_id
    hgnc.complete.set <- data.table::fread(file.path(data_dir, hgnc_file)) %>% as.data.frame()
    common_genes <- intersect(colnames(RTPT_tumor_mat), colnames(RTPT_CL_mat))
    common_genes <- intersect(common_genes, hgnc.complete.set$symbol)
    
    RTPT_tumor_mat <- RTPT_tumor_mat[,common_genes]
    RTPT_CL_mat <- RTPT_CL_mat[,common_genes]
    
    
    # 先处理hgnc.complete.set, 为什么要处理hgnc.complete.set呢？
    hgnc.complete.set <- filter(hgnc.complete.set, symbol %in% common_genes)
    hgnc.complete.set <- hgnc.complete.set[!duplicated(hgnc.complete.set$symbol),]
    rownames(hgnc.complete.set) <- hgnc.complete.set$symbol
    
    # 排序
    hgnc.complete.set <- hgnc.complete.set[common_genes,]
    
    # 转换成为ensembl id
    # 直接转换时有空值，所以必须先求交集，另外像AL627309.6 这样的基因其实是pseudogene。
    # colnames(RTPT_tumor_mat) <- hgnc.complete.set$ensembl_gene_id[match(colnames(RTPT_tumor_mat), hgnc.complete.set$symbol)]
    # colnames(RTPT_CL_mat) <- hgnc.complete.set$ensembl_gene_id[match(colnames(RTPT_CL_mat), hgnc.complete.set$symbol)]
    
    # 转换之前，gene id没有重复，交集为37035，转换之后 ensembl gene id 有重复。
    length(intersect(colnames(RTPT_tumor_mat),colnames(RTPT_CL_mat))) # 36979
    colnames(RTPT_tumor_mat) <- hgnc.complete.set$ensembl_gene_id
    colnames(RTPT_CL_mat) <- hgnc.complete.set$ensembl_gene_id
    
    # 测试重复
    summary(!duplicated(colnames(RTPT_tumor_mat)))
    length(!duplicated(colnames(RTPT_CL_mat)))
    RTPT_tumor_mat[1:2,37000:37035]
    RTPT_CL_mat[,37000:37035]
    # hgnc.complete.set[37000:37035,c("symbol","ensembl_gene_id")]
    # 再根据功能去掉一些gene: length(func_genes)=19651
    hgnc.complete.set[1:5,1:9]
    func_genes <- dplyr::filter(hgnc.complete.set, !locus_group %in% c('non-coding RNA', 'pseudogene'))$ensembl_gene_id
    head(func_genes)
    length(func_genes) # 19651
    
    length(intersect(colnames(RTPT_tumor_mat),colnames(RTPT_CL_mat))) # 36979
    genes_used <- intersect(colnames(RTPT_tumor_mat),colnames(RTPT_CL_mat))
    length(genes_used) # 36979
    genes_used <- intersect(genes_used, func_genes)
    length(genes_used) # 19641
    head(genes_used)
    RTPT_tumor_mat[1:5,1:15]
    dim(RTPT_tumor_mat)
    RTPT_tumor_mat = RTPT_tumor_mat[,genes_used]
    RTPT_CL_mat = RTPT_CL_mat[,genes_used]
    
    
    ################# for annotation file:
    RTPT_tumor_ann <- data.table::fread(file.path(data_dir, RootPath_Tumor_info_file)) %>% as.data.frame()
    RTPT_CL_ann <- data.table::fread(file.path(data_dir, RootPath_CL_info_file)) %>% as.data.frame()
    all_ann = rbind(ann,RTPT_tumor_ann,RTPT_CL_ann)
    RTPT_tumor_ann = dplyr::filter(all_ann, type=="rootpath")
    RTPT_CL_ann = dplyr::filter(all_ann, type=="rootpath")  # bug？
    
    
    
    return(list(RTPT_tumor_mat = RTPT_tumor_mat, RTPT_tumor_ann = RTPT_tumor_ann, RTPT_CL_mat = RTPT_CL_mat, RTPT_CL_ann = RTPT_CL_ann))
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
  seu_obj %<>% Seurat::FindClusters(resolution = global$mod_clust_res) # fix reductions bug
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

calc_gene_stats <- function(TCGA_mat, CCLE_mat, hgnc_file="database/hgnc_complete_set_7.24.2018.txt") {
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
# my_Tumor_file = args[1]
# my_info_file = args[2]
# dat = load_data(data_dir, RootPath_Tumor_file=my_Tumor_file, RootPath_Tumor_info_file=my_info_file)
dat = load_data(data_dir)
############################## start load my data
my_dat = load_my_data(data_dir)

# 有了dat之后就可以画图了：
####### for fig1 
######################################
# dplyr::mutate 是创建新变量； rbind是按照row粘贴到一起。
comb_ann <- rbind(
  dat$TCGA_ann %>% dplyr::select(sampleID, lineage, subtype, `Primary/Metastasis`) %>%
    dplyr::mutate(type = 'tumor'),
  dat$CCLE_ann %>% dplyr::select(sampleID, lineage, subtype, `Primary/Metastasis`) %>%
    dplyr::mutate(type = 'CL')
)

uncorrected_combined_type_plot = plot_uncorrected_data(dat$CCLE_mat, dat$TCGA_mat,comb_ann)
ggsave("Fig1.uncorrected_repeating.png", plot= uncorrected_combined_type_plot)

# add rootpath data Fig1:
comb_ann <- rbind(
  dat$TCGA_ann %>% dplyr::select(sampleID, lineage, subtype, `Primary/Metastasis`) %>% dplyr::mutate(type = 'tumor'),
  dat$CCLE_ann %>% dplyr::select(sampleID, lineage, subtype, `Primary/Metastasis`) %>% dplyr::mutate(type = 'CL'),
  dat$RTPT_tumor_ann %>% dplyr::select(sampleID, lineage, subtype, `Primary/Metastasis`) %>% dplyr::mutate(type = 'rootpath'),
  dat$RTPT_CL_ann %>% dplyr::select(sampleID, lineage, subtype, `Primary/Metastasis`) %>% dplyr::mutate(type = 'rootpath')
)

uncorrected_combined_type_plot_my_data = plot_uncorrected_customized_data(dat$CCLE_mat, dat$TCGA_mat, comb_ann, dat$RTPT_tumor_mat, dat$RTPT_CL_mat)
ggsave("Fig1.uncorrected_combined_type_plot.with.my.tumor.png", plot=uncorrected_combined_type_plot_my_data)
##################################################################################
# Fig2:
# type 是写到meta.data中
TCGA_obj <- create_Seurat_object(dat$TCGA_mat, dat$TCGA_ann, type='tumor')
CCLE_obj <- create_Seurat_object(dat$CCLE_mat, dat$CCLE_ann, type='CL')
TCGA_obj <- cluster_data(TCGA_obj)
CCLE_obj <- cluster_data(CCLE_obj)

# TCGA_my_obj <- create_Seurat_object(rbind(dat$TCGA_mat, dat$RTPT_tumor_mat), rbind(dat$TCGA_ann, dat$RTPT_tumor_ann), type='tumor')
# TCGA_my_obj <- cluster_data(TCGA_my_obj)
# CCLE_my_obj <- create_Seurat_object(rbind(dat$CCLE_mat, dat$RTPT_CL_mat), rbind(dat$CCLE_ann, dat$RTPT_CL_ann), type='tumor') #### bug
# CCLE_my_obj <- cluster_data(CCLE_my_obj)

tumor_DE_genes <- find_differentially_expressed_genes(TCGA_obj)
# tumor_DE_genes_1 <- find_differentially_expressed_genes(TCGA_my_obj)
CL_DE_genes   <- find_differentially_expressed_genes(CCLE_obj)
# CL_DE_genes_1 <- find_differentially_expressed_genes(CCLE_my_obj)

gene_stats   <- calc_gene_stats(dat$TCGA_mat, dat$CCLE_mat)

# gene_stats_1 <- calc_gene_stats(rbind(dat$TCGA_mat,dat$RTPT_tumor_mat), rbind(dat$CCLE_mat,dat$RTPT_CL_mat))
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
# DE_genes_1 <- full_join(tumor_DE_genes_1, CL_DE_genes_1, by = 'Gene', suffix = c('_tumor', '_CL')) %>%
#     mutate(tumor_rank = dplyr::dense_rank(-gene_stat_tumor),
#       CL_rank = dplyr::dense_rank(-gene_stat_CL),
#       best_rank = pmin(tumor_rank, CL_rank, na.rm=T)) %>%
#     dplyr::left_join(gene_stats_1, by = 'Gene')
# DE_gene_set_1 <- DE_genes_1 %>% dplyr::filter(best_rank < global$top_DE_genes_per) %>% .[['Gene']]
# cov_diff_eig_1 <- run_cPCA(TCGA_my_obj, CCLE_my_obj, global$fast_cPCA)
# how to compare differences between DE_genes and DE_genes_1 ??????


if(is.null(global$fast_cPCA)) {
    cur_vecs <- cov_diff_eig$vectors[, global$remove_cPCA_dims, drop = FALSE]
    # cur_vecs_1 <- cov_diff_eig_1$vectors[, global$remove_cPCA_dims, drop = FALSE]
} else {
    cur_vecs <- cov_diff_eig$rotation[, global$remove_cPCA_dims, drop = FALSE]
    # cur_vecs_1 <- cov_diff_eig_1$rotation[, global$remove_cPCA_dims, drop = FALSE]
}
# how to compare differences between cur_vecs and cur_vecs_1 ??????


# TCGA_my_mat = rbind(dat$TCGA_mat, dat$RTPT_tumor_mat)
# CCLE_my_mat = rbind(dat$CCLE_mat, dat$RTPT_CL_mat)
rownames(cur_vecs) <- colnames(dat$TCGA_mat)
# rownames(cur_vecs_1) <- colnames(TCGA_my_mat)

TCGA_cor <- resid(lm(t(dat$TCGA_mat) ~ 0 + cur_vecs)) %>% t()
CCLE_cor <- resid(lm(t(dat$CCLE_mat) ~ 0 + cur_vecs)) %>% t()
# TCGA_my_cor <- resid(lm(t(TCGA_my_mat) ~ 0 + cur_vecs_1)) %>% t()
# CCLE_my_cor <- resid(lm(t(CCLE_my_mat) ~ 0 + cur_vecs_1)) %>% t()

mnn_res   <- run_MNN(CCLE_cor, TCGA_cor,        k1 = global$mnn_k_tumor, k2 = global$mnn_k_CL, ndist = global$mnn_ndist, subset_genes = DE_gene_set)
# mnn_res_1 <- run_MNN(CCLE_my_cor, TCGA_my_cor,  k1 = global$mnn_k_tumor, k2 = global$mnn_k_CL, ndist = global$mnn_ndist, subset_genes = DE_gene_set_1)
# detected tied distances to neighbors, see ?'BiocNeighbors-ties'  ##### bug??????

combined_mat    <- rbind(mnn_res$corrected, CCLE_cor)
# combined_my_mat <- rbind(mnn_res_1$corrected, CCLE_my_cor)
comb_obj <- create_Seurat_object(combined_mat, comb_ann)
comb_obj <- cluster_data(comb_obj)
# Warning: The following arguments are not used: reduction
# comb_my_obj <- create_Seurat_object(combined_my_mat, comb_ann)
comb_my_obj <- cluster_data(comb_my_obj)

# fig2：
alignment =  Seurat::Embeddings(comb_obj, reduction = 'umap') %>% as.data.frame() %>% set_colnames(c('UMAP_1', 'UMAP_2')) %>% 
    rownames_to_column(var = 'sampleID') %>% left_join(comb_ann, by = 'sampleID')
# alignment2 = Seurat::Embeddings(comb_my_obj, reduction = 'umap') %>% as.data.frame() %>% set_colnames(c('UMAP_1', 'UMAP_2')) %>%
#     rownames_to_column(var = 'sampleID') %>% left_join(comb_ann, by = 'sampleID')

source(here::here('src', 'Figure2.R'))
fig_plot = Celligner_alignment_plot(alignment)
ggsave("Fig2.raw.umap.png", plot = fig_plot)

# DimPlot(comb_obj, reduction = 'umap',group.by = "orig.ident",tissue_colors)


selected_data <- alignment2[which(alignment2$type=='rootpath'),]
selected_data <- selected_data[!duplicated(selected_data$sampleID),]
fig_plot2 = Celligner_alignment_customized_plot(alignment2, selected_data)
ggsave("Fig2.all_comb.umap.png", plot = fig_plot2)
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


### 为了调试UMAP图，临时代码
# ggplot2::ggplot(uncorrected_alignment, ggplot2::aes(UMAP_1, UMAP_2)) +
#   ggplot2::geom_point(data = filter(uncorrected_alignment, type=='tumor'), alpha=0.6, size=0.5, pch=21, color='white', aes(fill=type)) +
#   ggplot2::geom_point(data = filter(uncorrected_alignment, type=='CL'), alpha=0.6, size=0.6, pch=3, aes(color=type), stroke=0.5) +
#   ggplot2::scale_color_manual(values=c(CL="#F8766D")) +
#   ggplot2::scale_fill_manual(values=c(tumor="#00BFC4")) +
#   ggplot2::xlab('UMAP 1') + ggplot2::ylab("UMAP 2") +
#   ggplot2::theme_classic() +
#   ggplot2::theme(legend.position='bottom', text = ggplot2::element_text(size=8),
#                  axis.text = ggplot2::element_text(size=6),
#                  axis.title = ggplot2::element_text(size=8),
#                  axis.ticks.length = unit(4,"lines"),
#                  legend.margin =ggplot2::margin(0,0,0,0),
#                  legend.box.margin=ggplot2::margin(-10,-30,-10,-30),
#                  axis.line = ggplot2::element_line(size = .9, color = "red", arrow=arrow(),linetype = 2, lineend = seq_along(3), linewidth = 2)
#                  )




# figure 3a

Celligner_aligned_data = Seurat::GetAssayData(comb_obj) %>% as.matrix() %>% t()
Celligner_aligned_my_data = Seurat::GetAssayData(comb_my_obj) %>% as.matrix() %>% t()
tumor_CL_cor = calc_tumor_CL_cor(Celligner_aligned_data, comb_ann)
tumor_CL_my_cor = calc_tumor_CL_cor(Celligner_aligned_my_data, comb_ann)
write.table(t(tumor_CL_cor) %>% as.data.frame() %>% rownames_to_column("ccle_name"), 
    file="tumor_CL_cor.transit.csv", sep=",", row.names = FALSE, quote = FALSE)
write.table(t(tumor_CL_my_cor) %>% as.data.frame() %>% rownames_to_column("ccle_name"), 
    file="tumor_CL_my_cor.transit.csv", sep=",", row.names = FALSE, quote = FALSE)
# select Rootpath tumor and root path cell line.
rootpath_tumor_cor = tumor_CL_my_cor[rownames(RootPath_CL_file),rownames(RootPath_Tumor_file)]
write.table(t(rootpath_rumor_cor) %>% as.data.frame() %>% rownames_to_column("ccle_name"),
  file="tumor_CL_select.transit.csv", sep=",", row.names= FALSE, quote =FALSE)

# fig3
get_cell_line_tumor_class <- function(tumor_CL_cor, alignment) {
  cl_tumor_classes <- apply(tumor_CL_cor, 2, function(x) cell_line_tumor_class(x, tumor_CL_cor, alignment)) %>% as.character()
  names(cl_tumor_classes) <- colnames(tumor_CL_cor)
  return(cl_tumor_classes)
}
source(here::here('src', 'Figure3.R'))
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
#save.image("main.Rdata")
