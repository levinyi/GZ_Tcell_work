library(stringr)
library(shiny)
library(here)
library(optparse)
library(ggplot2)
library(patchwork)
library(Seurat)
source(here('src', 'PlotFunctions.R'))
source(here('src', 'SeuratWrappers.R'))

srt_obj <- readRDS("integrated_seurat.rds")

d <- FetchData(object=srt_obj, vars=c(c("UMAP_1", "UMAP_2")))

tpk_gene <- function(srt, gene, assay) {
  counts <- srt@assays[[assay]]@counts
  dat <- counts[gene, ]
  tots <- colSums(counts)
  tpk <- 1000 * dat / tots
  return(tpk)
}

genes <- c("CD4", "CD8A", "FOXP3", "GAPDH")
assays <- c("SCT", "SCT", "SCT", "SCT")
for (igene in seq_along(genes)) {
  gene_tpk <- tpk_gene(srt_obj, genes[igene], assays[igene])
  new_col <- rep(NA, each=dim(d)[1])
  names(new_col) <- row.names(d)
  for (bc in names(new_col)) {
    if (bc %in% names(gene_tpk)) {
      new_col[bc] <- gene_tpk[bc]
    }
  }
  d[[genes[igene]]] <- new_col
}

# for single-sample analysis
# cluster_meta <- "SCT_snn_res.1.6"
# for multi-sample analysis
cluster_meta <- "integrated_snn_res.0.8"
phenos <- list("ID:CD8"=c(0,1,4,5,7),
               "ID:CD4"=c(2,6,8),
               "ID:CD4_Treg"=c(3))
for (ipheno in seq_along(phenos)) {
  new_col <- sapply(row.names(d),
                    function(bc) {
                      srt_obj[[cluster_meta]][bc,1] %in% phenos[ipheno][[1]]
                    })
  d[[names(phenos)[ipheno]]] <- new_col
}

write.csv(d, "InfiniTIL_data.csv")

