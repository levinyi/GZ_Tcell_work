library(Seurat)
library(ggplot2)
library(patchwork)
library(tidyverse)
library(here)

# script_dir <- "$10X/InfiniTIL.pipeline"
script_dir <- "/cygene2/pipeline/10X/InfiniTIL"
source(paste(script_dir, "src", "SeuratWrappers.R", sep = "/"))

## Load data

data_dir <- here('./')
srt_obj <- readRDS(paste0(data_dir, "/integrated_seurat.rds"))
vdj_data <- read.csv(paste0(data_dir, "/merged_tcr_table.csv"))

tcr_table <- select(vdj_data, c(barcode, chain, v_gene, j_gene, cdr3s_aa, TCR_clonality, raw_clonotype_id))
v_genes <- aggregate(tcr_table$v_gene, by = list("barcode"=tcr_table$barcode), FUN = paste, collapse = ",")
j_genes <- aggregate(tcr_table$j_gene, by = list("barcode"=tcr_table$barcode), FUN = paste, collapse = ",")
cdr3s <- aggregate(tcr_table$cdr3s_aa, by = list("barcode"=tcr_table$barcode), FUN = unique)
counts <- aggregate(tcr_table$TCR_clonality, by = list("barcode"=tcr_table$barcode), FUN = unique)
orig <- aggregate(tcr_table$raw_clonotype_id, by = list("barcode"=tcr_table$barcode), FUN = unique)
new_table <- cbind(v_genes$barcode, v_genes$x, j_genes$x, cdr3s$x, counts$x, orig$x)
colnames(new_table) <- c("barcode", "v_genes", "j_genes", "cdr3s_aa", "TCR_clonality", "original_clone")
new_table <- as.data.frame(new_table)

### NeoTCR4 AuCell scoring

library(AUCell)
library(GSEABase)

score_aucell <- function(srt_obj, gene_list, drop_f=0.05) {
    raw_matrix <- srt_obj@assays$RNA@counts
    ranking <- AUCell_buildRankings(raw_matrix, plotStats = T, verbose = T, splitByBlocks=TRUE)
    max_rank <- ceiling(quantile(srt_obj$nFeature_RNA, drop_f))
    print(paste("Max rank:", max_rank))
    # build genesets
    module_name <- names(gene_list)[1]
    genesets <- list("set" = GeneSet(gene_list[[module_name]], setName = "set"))
    cells_auc <- AUCell_calcAUC(geneSets = genesets[["set"]],
                                rankings = ranking,
                                aucMaxRank = max_rank)
    srt_obj <- AddMetaData(srt_obj,
                           metadata = cells_auc@assays@data$AUC[1, ],
                           col.name = paste0(module_name, "_AUCscore"))
    assignments <-  AUCell_exploreThresholds(cells_auc, assignCells = T)
    return(srt_obj)
}

library(readxl)
full_table <- read_excel(paste0("/cygene2/pipeline/10X/InfiniTIL/index/science.abl5447_table_s4.xlsx"), sheet = "Signatures from this study")
Neo_CD4 <- full_table$NeoTCR4
Neo_CD4 <- Neo_CD4[!is.na(Neo_CD4)]
Neo_CD8 <- full_table$NeoTCR8
Neo_CD8 <- Neo_CD8[!is.na(Neo_CD8)]

neocd4_name <- "NeoTCR4"
gene_list <- list(Neo_CD4)
names(gene_list) <- neocd4_name
srt_obj <- score_aucell(srt_obj, gene_list)

print(FeaturePlot(srt_obj, features = c(paste0(neocd4_name, "_AUCscore"))))

new_col <- rep(NA, each=dim(new_table)[1])
names(new_col) <- new_table$barcode
for (bc in names(new_col)) {
  if (bc %in% colnames(srt_obj)) {
    new_col[bc] <- srt_obj[[paste0(neocd4_name, "_AUCscore")]][bc, 1]
  }
}
new_table[[neocd4_name]] <- new_col

neocd8_name <- "NeoTCR8"
gene_list <- list(Neo_CD8)
names(gene_list) <- neocd8_name
srt_obj <- score_aucell(srt_obj, gene_list)

print(FeaturePlot(srt_obj, features = c(paste0(neocd8_name, "_AUCscore"))))

new_col <- rep(NA, each=dim(new_table)[1])
names(new_col) <- new_table$barcode
for (bc in names(new_col)) {
  if (bc %in% colnames(srt_obj)) {
    new_col[bc] <- srt_obj[[paste0(neocd8_name, "_AUCscore")]][bc, 1]
  }
}
new_table[[neocd8_name]] <- new_col

write_csv(new_table, paste0(data_dir, "/TCR_stats_barcodes.csv"))
