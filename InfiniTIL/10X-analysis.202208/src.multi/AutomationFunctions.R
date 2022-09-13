# Helper and wrapper functions for automated analysis of 10X data.
# This R file is meant to be sourced in the context of a call to BatchAnalysis

# Author: Seon Kinrot (RootPath Genomics)

# dependancies
library(stringr)
library(shiny)
library(here)
library(ggplot2)
library(patchwork)
library(Seurat)
source(here('src', 'PlotFunctions.R'))
source(here('src', 'SeuratWrappers.R'))

# global list of gene lists/ modules
GLOBAL_MODULES <- list(
  "YostEx" = c("CTSW", "GNLY", "PRF1", "SLA2", "HAVCR2", "GZMB", "ITGAE",
               "GALNT2", "ACP5", "CXCR6", "ENTPD1", "KRT86", "TIGIT", "LAYN",
               "JAML", "AC092580.4", "AHI1", "ALOX5AP", "SYNGR2", "CARD16",
               "FKBP1A", "GOLIM4", "SNX9", "TNFRSF18", "ZBED2", "TNFRSF9",
               "VCAM1", "CXCL13", "GEM"),
  "TCRx_sig" = c("CXCL13", "CTLA4", "DUSP4", "ENTPD1", "RBPJ", "TNFSF4", "GK",
                 "PDCD1", "LAG3", "FKBP1A", "MYO7A", "FABP5", "GAPDH", "TIGIT",
                 "CDCA7", "ZBED2", "RGS1", "HAVCR2", "MKI67", "GZMB", "SPRY1",
                 "PHLDA1", "ASPM", "MYO1E", "TOX", "LSP1", "MCM10", "MAD2L1",
                 "OGN", "PCLAF", "UBE2T", "GNG4", "LY6E", "TPI1", "TNFRSF18",
                 "RGS2", "CLEC2D", "SAMSN1", "TSHZ2", "HMGN2"),
  "TR_1.0" = c("CXCL13", "TNFRSF9", "HAVCR2", "CTLA4", "ENTPD1",
               "TIGIT", "VCAM1", "TNFRSF18", "PHLDA1", "LAYN", "DUSP4",
               "KRT86", "GEM", "CD82", "RBPJ", "PDCD1", "MIR155HG",
               "SNX9", "LAG3", "TOX", "CD27", "ACP5", "NDFIP2", "CCL3",
               "SIRPG", "RGS2", "LYST", "FAM3C", "IGFLR1", "GZMB"),
  "NeoTCR4" = c("CXCL13", "HMOX1", "ETV7", "ADGRG1", "PDCD1", "ENTPD1",
                "CCDC50", "TOX", "CD4", "TIGIT", "TNFRSF18", "NMB", "MYL6B",
                "AHI1", "MAF", "IFNG", "LAG3", "CXCR6", "IGFLR1", "DUSP4",
                "ACP5", "LINC01943", "LIMS1", "BATF", "PCED1B", "ITGAL",
                "YPEL2", "MAL", "PPT1", "ELMO1", "MIS18BP1", "TMEM173", "ADI1",
                "SLA", "GALM", "LBH", "SECISBP2L", "CTSB", "C17orf49",
                "CORO1B"),
  "NeoTCR8" = c("ATP10D", "GZMB", "ENTPD1", "KIR2DL4", "LAYN", "HTRA1", "CD70",
                "CXCR6", "HMOX1", "ADGRG1", "LRRN3", "ACP5", "CTSW", "GALNT2",
                "LINC01480", "CARS", "LAG3", "TOX", "PTPRCAP", "ASB2", "ITGB7",
                "PTMS", "CD8A", "GPR68", "NSMCE1", "ABI3", "SLC1A4", "PLEKHF1",
                "CD8B", "LINC01871", "CCL4", "NKG7", "CLIC3", "NDFIP2",
                "PLPP1", "PCED1B", "CXCL13", "PDCD1", "PRF1", "HLA-DMA",
                "GPR25", "CD9", "TIGIT", "HLA-DRB5", "SYTL3", "SLF1", "NEK1",
                "CASP1", "SMC4", "TSEN54", "PLSCR1", "GNPTAB", "HLA-DPB1",
                "PLEKHA1", "ARHGAP9", "ALOX5AP", "SH3BP1", "NCF4", "NELL2",
                "GATA3", "PPM1M", "TNFRSF1A", "AC022706.1", "MCM5", "HLA-DRB1",
                "TNFSF10", "TRIM21", "HDLBP", "ERN1", "CALHM2", "SASH3",
                "ACTA2", "MAST4", "CAPG", "MPST", "IGFLR1", "GZMA", "CD27",
                "ITGAE", "SLA2", "RHOC", "COMMD8", "MYO1G", "SP140", "PHPT1",
                "CD2BP2", "PLEKHO1", "STAM", "MRPL16", "IL2RB", "ID2",
                "TESPA1", "GOLGA8B", "MIS18BP1", "VAMP5", "DAPK2", "HLA-DPA1",
                "TSG101", "IL4R", "CCND2", "CTSC", "TRAF3IP3", "NLRC3",
                "ORAI3", "GNLY", "MIR155HG", "CARD16", "CD82", "ECH1", "JAML",
                "EEF1G", "ETFB", "DAXX", "RBM4", "HCST", "RAB27A", "YPEL2",
                "CHST12", "ARPC1B", "PDIA4", "PDIA6", "AC243960.1", "TBC1D10C",
                "PTPN6", "PYCARD", "BST2", "BTN3A2", "MTG1", "MLEC", "DUSP4",
                "GSDMD", "SLAMF1", "IFI6", "PCID2", "GIMAP1", "ITGA1",
                "CSNK2B", "CDK2AP2", "MYO1F", "AC004687.1", "PTTG1",
                "APOBEC3C", "TSPAN14", "MOB3A", "STXBP2", "LCP2", "PLA2G16",
                "LINC00649", "CST7", "TADA3", "SIT1", "APOBEC3G", "SUSD3",
                "CD3G", "CCL5", "CDC25B", "TNFRSF1B", "HMGN3", "THEMIS",
                "ASF1A", "CTNNB1", "FIBP", "CCDC85B", "POLR3GL", "GIMAP6",
                "ARL6IP1", "CALCOCO2", "CCPG1", "KLRB1", "ACAA2", "ISG15",
                "EIF4A1", "CAT", "MANF", "XAB2", "GRINA", "GLO1", "LSM2",
                "SLFN5", "FKBP1A", "AKNA", "TAP1", "LMO4", "APEH", "C12orf75",
                "TMEM14A", "DNPH1", "C17orf49", "NUDT5", "MGAT1", "CCDC69",
                "EIF4EBP1", "PDHB", "ARL3", "UCP2", "IFI35", "HSBP1", "LYST",
                "MRFAP1L1", "ITGAL", "AIP", "RASAL3", "CAPN1", "ITGB1", "RBPJ",
                "LBH", "DYNLL1", "NME2", "MT1F", "SYNGR2", "ABTB1", "ZGPAT",
                "CD63", "ILK", "SKA2", "TMEM204", "ACO2", "HOPX", "CRIP1",
                "OXNAD1", "CCS", "GRAP2", "GSTO1", "HADHB", "IL16", "PIN4",
                "CUEDC2", "CALM3", "SAMSN1", "HM13", "SNAP23", "LPCAT4",
                "FAAP20", "EFHD2", "PRDX3", "CCM2", "C22orf39", "SDHA",
                "ARRDC1", "MAP4K1", "NDUFA13", "IL27RA", "C14orf119")
  )

#' Function to generate a string, describing how much time
#' has past since start_time
time_string <- function(start_time) {
  now <- Sys.time()
  elapsed <- as.numeric(now) - as.numeric(start_time)
  secs <- as.integer(elapsed) %% 60
  mins <- as.integer((elapsed %% 3600) / 60)
  hrs <- as.integer(elapsed / 3600)
  if (hrs == 0) {
    string <- paste0(mins, ":", secs, " minutes")
  } else {
    string <- paste0(hrs, "::", mins, ":", secs, " hours")
  }
  return(string)
}

#' Function to handle and check CLI input arguments for
#' validity and consistency. Returns a list of verified inputs,
#' or prints an error message and returns NULL
handle_options <- function(parsed_options, args=NULL) {
  in_list <- list()

  # start by treating evry argument as "simple" - if it isn't, modify it later
  # handle inputs - simple args
  for (arg in args) {
    in_list[[arg]] <- parsed_options[[arg]]
  }

  # handle inputs - data directory
  if ("input" %in% args) {
    data_dir <- parsed_options$input
    if (!dir.exists(data_dir)) {
      cat("Error -", data_dir, "does not exist!\n", sep = " ")
      return(NULL)
    }
    in_list$data_dir <- data_dir
  }
  # handle inputs - save directory
  if ("output" %in% args) {
    save_dir <- parsed_options$output
    if (is.null(save_dir)) {
      save_dir <- paste(data_dir, "analysis", sep="/")
    }
    if (!dir.exists(save_dir)) {
      dir.create(save_dir)
    }
    in_list$save_dir <- save_dir
  }

  # allow "shortcut" so user doesn't have to specify both -c and -s
  if ("has_csp" %in% args & "sep_csp" %in% args) {
    if ((length(in_list$sep_csp) > 0) & in_list$sep_csp) {
      in_list$has_csp <- TRUE
    }
  }
  # handle inputs - species
  if ("mouse" %in% args) {
    species <- "human"
    if ((length(parsed_options$mouse) > 0) & parsed_options$mouse) {
      species <- "mouse"
    }
    in_list$species <- species
  }

  return(in_list)
}

#' Find all GEX R1 fastq files in a folder
gex_r1_list <- function(directory) {
  gex_pattern <- "GEX.*_R1_.*fastq[gz.]*$"
  fq_files <- list.files(in_list$data_dir, pattern = gex_pattern)
  return(fq_files)
}

#' wrapper for file names
handle_pos_args <- function(parsed_args, in_list) {
  in_list$project_names <- parsed_args
  # no project names provided - assign numerically to each detected project
  if (length(in_list$project_names)==0) {
    # detect by FASTQ files
    fq_files <- gex_r1_list(in_list$data_dir)
    if (length(fq_files) > 0) {
      in_list$project_names <- sapply(X = fq_files,
                                      FUN = function(fl) {
                                        gex_name <- strsplit(fl, "_S1_")[[1]][1]
                                        project_name <- sub("[-_]*GEX[-_]*", "",
                                                            gex_name)
                                        return(project_name)
                                      }
                                    )
    }
    else {  # detect by directories
      subdirs <- list.dirs(in_list$data_dir,
                           recursive = FALSE,
                           full.names = FALSE)
      project_dirs <- grep("_GEX", subdirs, value = TRUE)
      in_list$project_names <- sapply(X = project_dirs,
                                      FUN = function(x) {
                                        return(sub("_GEX*", "", x))
                                      }
                                    )
    }
  }  # end automated project name inference
  return(in_list)
}

#' Make sure everything is ready to run CellRanger. Returns True/ False
#' as an "answer" to this question
validate_cr <- function(in_list) {
  # CellRanger not installed or not in path
  status <- invisible(system("which cellranger"))
  if (status != 0) {
    cat("Error - cellranger not installed or not in system path")
    return(FALSE)
  }
  # no transcriptome file
  if (!dir.exists(in_list$transcriptome)) {
    cat("Error -", transcriptome, "does not exist!\n", sep = " ")
    return(FALSE)
  }
  # no vdj_ref file
  if (in_list$has_vdj & !dir.exists(in_list$vdj_ref)) {
    cat("Error -", vdj_ref, "does not exist!\n", sep = " ")
    return(FALSE)
  }
  # no feature_ref file
  if (in_list$has_csp & !file.exists(in_list$feature_ref)) {
    cat("Error -", feature_ref, "does not exist!\n", sep = " ")
    return(FALSE)
  }

  return(TRUE)
}

#' Automatically generate sample name from fastq file name
gex_sample_name <- function(directory, iproject) {
  fq_file <- gex_r1_list(directory)[iproject]
  gex_sample <- strsplit(fq_file, "_S1_")[[1]][1]
  return(gex_sample)
}

#' Function for generating a library file for a sample handled by CR
write_lib_file <- function(in_list, iproject) {
  lib_file <- paste0(in_list$data_dir, "/libs", iproject, ".csv")
  fastqs <- c(in_list$data_dir)
  gex_sample <- gex_sample_name(in_list$data_dir, iproject)
  samples <- c(gex_sample)
  types <- "Gene Expression"
  # add CSP data
  if (in_list$has_csp) {
    fastqs <- c(fastqs, in_list$data_dir)
    csp_sample <- sub("GEX", "CSP", gex_sample)
    samples <- c(samples, csp_sample)
    types <- c(types, "Antibody Capture")
  }
  file_out <- as.data.frame(list("fastqs"=fastqs,
                                 "sample"=samples,
                                 "library_type"=types))
  write.csv(file_out, file = lib_file,
            row.names = FALSE, col.names = TRUE)
  return(lib_file)
}

#' Wrapper function for running CellRanger on a single sample
run_cr_sample <- function(in_list, iproject) {
  # generate library file
  lib_file <- write_lib_file(in_list, iproject)

  # id for CellRanger run
  cr_id <- "GEX"
  if (in_list$has_csp) {
    cr_id <- "GEX-CSP"
  }
  gex_sample <- gex_sample_name(in_list$data_dir, iproject)
  cr_id <- sub("GEX", cr_id, gex_sample)

  # run CellRanger count
  count_cmd <- paste("cellranger", "count", "--id", cr_id,
                     "--transcriptome", in_list$transcriptome,
                     "--feature-ref", in_list$feature_ref,
                     "--libraries", lib_file,
                     "--localcores", in_list$cpus,
                     "--localmem",in_list$max_memory)
  system(count_cmd)

  # CellRanger vdj portion
  if (in_list$has_vdj) {
    vdj_sample <- sub("GEX", "VDJ", gex_sample)
    sample <- c(vdj_sample)
    vdj_cmd <- paste("cellranger", "vdj", "--id", sample,
                     "--reference", in_list$vdj_ref,
                     "--fastqs", in_list$data_dir,
                     "--sample", sample,
                     "--localcores", in_list$cpus,
                     "--localmem",in_list$max_memory)
    system(vdj_cmd)
  }

  return(NULL)
}

#' Wrapper function to load all post-CR data to Seurat
load_all_data <- function(in_list) {
  n_projects <- length(in_list$project_names)
  gex_list <- vector(mode = "list", length = n_projects)
  vdj_list <- vector(mode = "list", length = n_projects)
  csp_list <- vector(mode = "list", length = n_projects)

  for (iproject in seq_len(n_projects)) {
    gex_dir <- paste0(in_list$data_dir, "/",
                      in_list$project_names[iproject], "_GEX")
    if (in_list$has_csp & !in_list$sep_csp) { # GEX and CSP are fused
      gex_dir <- sub("GEX", "GEX-CSP", gex_dir)
    }
    gex_list[[iproject]] <- suppressWarnings(
                              load_seurat(
                                prj_dir = gex_dir,
                                prj_nm = in_list$project_names[[iproject]])
                              )

    if (in_list$has_vdj) {
      vdj_dir <- paste0(in_list$data_dir, "/",
                        in_list$project_names[iproject], "_VDJ")
      vdj_list[[iproject]] <- load_vdj(tcr_dir = vdj_dir)
    }

    if (in_list$has_csp) {
      if (in_list$sep_csp) {
        # load CSP data
        csp_dir <- sub("GEX", "CSP", gex_dir)
        csp_list[[iproject]] <- suppressWarnings(
                                  load_seurat(prj_dir = csp_dir,
                                    prj_nm = in_list$project_names[[iproject]])
                                  )
        # add to GEX. It shouldn't matter whether CSP was separate any more
        gex_list[[iproject]] <- add_assay(main_srt = gex_list[[iproject]],
                                          secondary_srt = csp_list[[iproject]],
                                          assay = "RNA",
                                          new_name = "ADT",
                                          filter_main = TRUE)
      }
      # CSP data is loaded in GEX
      csp_obj <- gex_list[[iproject]]
      # remove assays so they can be fused back later with the original names
      gex_list[[iproject]][["ADT"]] <- NULL
      DefaultAssay(csp_obj) <- "ADT"
      csp_obj[["RNA"]] <- NULL
      csp_list[[iproject]] <- csp_obj
    }  # end has_csp
  }  # end project loop

  # set names
  names(gex_list) <- in_list$project_names
  names(csp_list) <- in_list$project_names
  names(vdj_list) <- in_list$project_names

  data_list <- list("gex_list"=gex_list,
                    "csp_list"=csp_list,
                    "vdj_list"=vdj_list)
  return(data_list)
}

#' Wrapper function for running QC through Shiny and returning
#' the post-QC list of GEX objects
run_qc_shiny <- function(host=host, gex_list, shiny_dir, in_list=NULL, port=NA) {
  if (is.na(port)) {
    port <- getOption("shiny.port")
  }
  qc_dir <- paste0(shiny_dir, "/data")
  saveRDS(gex_list, paste0(qc_dir, "/input.rds"))
  runApp(host=host, shiny_dir, launch.browser = TRUE, port = port)
  stopApp() # shouldn't be necessary, but just in case
  if (!file.exists(paste0(qc_dir, "/output.rds"))) {
    cat("Failed to detect post-QC data; reverting to auto-QC\n")
    if (!is.null(in_list)) {
      return(run_auto_qc(gex_list, in_list))
    }
    return(NULL)
  } else {
    gex_list <- readRDS(paste0(qc_dir, "/output.rds"))
    return(gex_list)
  }
}

#' Wrapper for coarse, automatic QC
run_auto_qc <- function(gex_list, in_list) {
  n_projects <- length(gex_list)
  for (iproject in seq_len(n_projects)) {
    qc_res <- run_qc(srt_obj=gex_list[[iproject]],
                     #fsets are %mito and %ribo
                     min_fsets=c(0,0), max_fsets=c(10, 100),
                     min_feats=200, max_feats=Inf,
                     min_umi=0, max_umi=Inf,
                     species=in_list$species,
                     return_plot=TRUE, return_object=TRUE)
    gex_list[[iproject]] <- qc_res$filtered_seurat
    plot_list <- qc_res$plots
    qc_plot <- plot_list[[1]]
    for (panel in qc_plot[2:length(qc_plot)]) {
      qc_plot <- qc_plot + panel
    }
    ggsave(paste0(save_dir, "/", project_names[[iproject]], "_autoQC.png"),
           plot = qc_plot ,width = 8, height = 6)
  }
  return(gex_list)
}

#' Wrapper for adding CSP assays to GEX objects
#' while retaining as much information as possible
fuse_cite_seq <- function(integrated_gex, integrated_csp, integrated=TRUE) {
  integrated_gex <- add_assay(main_srt = integrated_gex,
                              secondary_srt = integrated_csp,
                              assay = "ADT",
                              new_name = "ADT",
                              add_reductions = TRUE,
                              filter_main = TRUE)
  csp_assay <- "ADT"
  gex_assay <- "RNA"
  if (integrated) {
    integrated_gex <- add_assay(main_srt = integrated_gex,
                                secondary_srt = integrated_csp,
                                assay = "integrated",
                                new_name = "adt_integrated",
                                add_reductions = TRUE,
                                filter_main = TRUE)
    csp_assay <- "adt_integrated"
    gex_assay <- "integrated"
  }
  # cluster CSP data
  DefaultAssay(integrated_gex) <- csp_assay
  # without integration, need to first normalize the data
  if (!integrated) {
    integrated_gex <- normalize_and_find_var_feats(
      integrated_gex, normalization_method = "CLR", do_scale = FALSE)
  }
  integrated_gex <- ScaleData(object = integrated_gex, verbose = FALSE,
                              vars.to.regress = "nCount_ADT")
  integrated_gex <- suppressWarnings(cluster_seurat(srt_obj = integrated_gex))
  # prevent GEX dimensionality reductions from overwriting the CSP ones later
  for (red in c("pca", "umap")) {
    # key change to ADTUMAP_ and ADTPC_
    integrated_gex[[red]]@key <- paste0("ADT", integrated_gex[[red]]@key)
    colnames(integrated_gex[[red]]@cell.embeddings) <-
      paste0("ADT", colnames(integrated_gex[[red]]))
    # copy reduction and add "adt_" prefix (e.g. "umap" becomes "adt_umap")
    temp <- integrated_gex[[red]]
    integrated_gex[[red]] <- NULL
    integrated_gex[[paste0("adt_", red)]] <- temp
  }
  rm(temp)
  DefaultAssay(integrated_gex) <- gex_assay

  return(integrated_gex)
}

#' Wrapper function for generating all "standard" plots
batch_plots <- function(integrated_gex, in_list, gene_modules=NULL) {
  n_projects <- length(in_list$project_names)
  # separate UMAP for each sample of origin
  samples_gex_plot <- DimPlot(object = integrated_gex, split.by = "orig.ident",
                              group.by = "orig.ident") +
    ggtitle("GEX UMAP per sample of origin")
  ggsave(paste0(in_list$save_dir, "/GEX_UMAP_by_sample.png"),
         samples_gex_plot,
         width = 4*n_projects, height = 6)
  # separate UMAP for cells predicted to belong to each cell cycle phase
  temp_var <- tryCatch(
    {cc_gex_plot <- DimPlot(object = integrated_gex, split.by = "Phase",
                            group.by = "Phase") +
      ggtitle("GEX UMAP per cell cycle phase")
    ggsave(paste0(in_list$save_dir, "/GEX_UMAP_by_ccPhase.png"),
           cc_gex_plot,
           width = 16, height = 6)
    },
    error = function(cond) {return(NULL)}
  )
  # graph-based clusters at increasing resolution
  integrated <- n_projects > 1
  gex_prefix <- "integrated"
  csp_prefix <- "adt_integrated"
  if (!integrated) {
    gex_prefix <- "SCT"
    csp_prefix <- "ADT"
  }
  gex_cluster_plt <- DimPlot(object = integrated_gex,
                             group.by = paste0(gex_prefix, "_snn_res.0.4"),
                             label = T, label.box = T, pt.size = 1) +
    DimPlot(object = integrated_gex,
            group.by = paste0(gex_prefix, "_snn_res.0.8"),
            label = T, label.box = T, pt.size = 1) +
    DimPlot(object = integrated_gex,
            group.by = paste0(gex_prefix, "_snn_res.1.2"),
            label = T, label.box = T, pt.size = 1)
  ggsave(filename = paste0(in_list$save_dir, "/GEX_UMAP_GexClusterByRes.png"),
         plot = gex_cluster_plt,
         width = 20,
         height = 6)
  # plot clustering for CSP data
  if (in_list$has_csp) {
    csp_cluster_plt <- DimPlot(object = integrated_gex,
                               group.by = paste0(csp_prefix, "_snn_res.0.4"),
                               reduction = "adt_umap",
                               label = T, label.box = T, pt.size = 1) +
      DimPlot(object = integrated_gex,
              group.by = paste0(csp_prefix, "_snn_res.0.8"),
              reduction = "adt_umap",
              label = T, label.box = T, pt.size = 1) +
      DimPlot(object = integrated_gex,
              group.by = paste0(csp_prefix, "_snn_res.1.2"),
              reduction = "adt_umap",
              label = T, label.box = T, pt.size = 1) +
      DimPlot(object = integrated_gex,
              group.by = paste0(gex_prefix, "_snn_res.1.6"),
              label = T, label.box = T, pt.size = 1) +
      DimPlot(object = integrated_gex,
              group.by = paste0(gex_prefix, "_snn_res.2"),
	      label = T, label.box = T, pt.size = 1)
    ggsave(filename = paste0(in_list$save_dir, "/CSP_UMAP_CspClusterByRes.png"),
           plot = csp_cluster_plt,
           width = 18,
           height = 12)
  }
  # plot TCR clonality on GEX UMAP
  if (in_list$has_vdj) {
    path <- paste0(in_list$save_dir, "/GEX_UMAP_TCRclonality.png")
    gex_clone_plot <- plot_clonality(integrated_gex,
                                     save_path = path)
  }
  # feature plots for specific genes
  if (in_list$species == "human") {
    # add modules to score for potential TCR reactivity prediction
    gene_modules <- c(GLOBAL_MODULES, gene_modules)
    module_names <- paste0(names(gene_modules), "_AUCscore")
    integrated_gex <- score_modules(srt_obj = integrated_gex,
                                    gene_list = gene_modules,
                                    module_name = module_names)
    plot_mods <- c(rootpath_jet_scaling(value_range = c(0, 0.25)))
    module_plot_panel <- multi_feature_plot(srt_obj = integrated_gex,
                                            feature_list = module_names,
                                            ggplot_mods = plot_mods)
    ggsave(filename = paste0(in_list$save_dir,
                             "/GEX_UMAP_ModuleScoresAUC.png"),
           plot = module_plot_panel,
           width = 12, height = 8)
    # markers from Yost et al.
    yost_markers <- c("CD3D","CD8A", "CD4", "FOXP3", "CCR7", "IL26", "CD200",
                     #"EOMES", "KLRD1", "IFNG", "HAVCR2")
		      "EOMES", "KLRD1", "IFNG", "HAVCR2", "CXCL13", "ENTPD1", "PDCD1", "ITGAE")
    yost_markers <- yost_markers[yost_markers %in% rownames(integrated_gex)]
    yost_fp <- multi_feature_plot(srt_obj = integrated_gex,
                                  feature_list = yost_markers,
                                  ggplot_mods=c(rootpath_jet_scaling())
                                )
    ggsave(filename = paste0(in_list$save_dir, "/GEX_UMAP_YostFP.png"),
           plot = yost_fp,
           width = 18, height = 12)
    # set of genes identified from literature to be plotted individually
    plot_dir <- paste0(in_list$save_dir, "/Feature_plots_GEX/")
    if (!dir.exists(plot_dir)) {
      dir.create(plot_dir)
    }
    suppressWarnings(gen_feature_plots(gex_seurat = integrated_gex,
                                       plot_dir = plot_dir))
  }

  return(integrated_gex)
}

#' Wrapper function to add columns to TCR tables
tcr_selection_table <- function(vdj_data, srt_obj, features=NULL, assays=NULL,
                                metadata=NULL, save_path=NULL) {
  # construct base table
  library(tidyverse)
  tcr_table <- dplyr::select(vdj_data,
                             c(barcode, chain, v_gene, j_gene, cdr3s_aa,
                               TCR_clonality, original_clone_id)
                             )
  v_genes <- aggregate(tcr_table$v_gene,
                       by = list("barcode"=tcr_table$barcode),
                       FUN = paste, collapse = ",")
  j_genes <- aggregate(tcr_table$j_gene,
                       by = list("barcode"=tcr_table$barcode),
                       FUN = paste, collapse = ",")
  cdr3s <- aggregate(tcr_table$cdr3s_aa,
                     by = list("barcode"=tcr_table$barcode),
                     FUN = unique)
  counts <- aggregate(tcr_table$TCR_clonality,
                      by = list("barcode"=tcr_table$barcode),
                      FUN = unique)
  orig <- aggregate(tcr_table$original_clone_id,
                    by = list("barcode"=tcr_table$barcode),
                    FUN = unique)
  new_table <- cbind(v_genes$barcode, v_genes$x, j_genes$x,
                     cdr3s$x, counts$x, orig$x)
  colnames(new_table) <- c("barcode", "v_genes", "j_genes",
                           "cdr3s_aa", "TCR_clonality", "original_clone")
  new_table <- as.data.frame(new_table)

  # add TPK for listed features
  if (length(features) != length(assays)) {
    cat("Features and assays do not match in length;",
        "skipping feature-based columns\n")
  }
  else {
    for (ifeat in seq_along(features)) {
      gene_tpk <- tpk_gene(srt_obj, features[ifeat], assays[ifeat])
      new_col <- rep(NA, each=dim(new_table)[1])
      names(new_col) <- new_table$barcode
      for (bc in names(new_col)) {
        if (bc %in% names(gene_tpk)) {
          new_col[bc] <- gene_tpk[bc]
        }
      }
      col_name <- paste0(features[ifeat], "_TPK")
      new_table[[col_name]] <- new_col
    }
  }

  # add specific metadata
  for (imeta in seq_along(metadata)) {
    new_col <- rep(NA, each=dim(new_table)[1])
    names(new_col) <- new_table$barcode
    for (bc in names(new_col)) {
      if (bc %in% colnames(srt_obj)) {
        new_col[bc] <- srt_obj[[metadata[imeta]]][bc, 1]
      }
    }
    new_table[[metadata[imeta]]] <- new_col
  }

  if (is.null(save_path)) {
    return(new_table)
  }
  write.csv(new_table, file = save_path)
  return(save_path)
}
