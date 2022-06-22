# Helper and wrapper functions for automated analysis of 10X data.
# This R file is meant to be sourced in the context of a call to BatchAnalysis

# Author: Seon Kinrot (RootPath Genomics)

# dependancies
library(stringr)
library(shiny)
library(here)
library(optparse)
library(ggplot2)
library(patchwork)
library(Seurat)
source(here('src', 'PlotFunctions.R'))
source(here('src', 'SeuratWrappers.R'))
# **Note**
# The following libraries are dependancies, but should already be loaded if
# this file is called from BatchAnalysis:
# SeuratWrappers (which also loads Seurat, ggplot2 and patchwork)
# optparse
# here (starting at the home directory of this repository)
# ---------
# The code was modified to include those dependancies explicity

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
    fq_files <- gex_r1_list(in_list$data_dir)
    in_list$project_names <- sapply(X = fq_files,
                                    FUN = function(fl) {
                                      gex_name <- strsplit(fl, "_S1_")[[1]][1]
                                      project_name <- sub("[-_]*GEX[-_]*", "",
                                                          gex_name)
                                      return(project_name)
                                    })
  }
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
run_qc_shiny <- function(host,gex_list, shiny_dir, in_list=NULL, port=NA) {
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
    plot_list <- qc_res$plot_list
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
batch_plots <- function(integrated_gex, in_list) {
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
            label = T, label.box = T, pt.size = 1) +
    DimPlot(object = integrated_gex,
            group.by = paste0(gex_prefix, "_snn_res.1.6"),
            label = T, label.box = T, pt.size = 1) +
    DimPlot(object = integrated_gex,
            group.by = paste0(gex_prefix, "_snn_res.2"),
            label = T, label.box = T, pt.size = 1)
  ggsave(filename = paste0(in_list$save_dir, "/GEX_UMAP_GexClusterByRes.png"),
         plot = gex_cluster_plt,
         width = 20,
         height = 12)
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
              label = T, label.box = T, pt.size = 1)
    ggsave(filename = paste0(in_list$save_dir, "/CSP_UMAP_CspClusterByRes.png"),
           plot = csp_cluster_plt,
           width = 20,
           height = 6)
  }
  # plot TCR clonality on GEX UMAP
  if (in_list$has_vdj) {
    path <- paste0(in_list$save_dir, "/GEX_UMAP_TCRclonality.png")
    gex_clone_plot <- plot_clonality(integrated_gex,
                                     save_path = path)
  }
  # feature plots for specific genes
  if (in_list$species == "human") {
    # T-cell exhaustion score from Yost et al.
    # add Yost exhaustion score
    if (integrated) {
      var_feats <- integrated_gex@assays$integrated@var.features
    } else {
      var_feats <- integrated_gex@assays$SCT@var.features
    }
    yost_exhaustion <- c("CTSW", "GNLY", "PRF1", "SLA2", "HAVCR2", "GZMB",
                         "ITGAE", "GALNT2", "ACP5", "CXCR6", "ENTPD1", "KRT86",
                         "TIGIT", "LAYN", "JAML", "AC092580.4", "AHI1",
                         "ALOX5AP", "SYNGR2", "CARD16", "FKBP1A", "GOLIM4",
                         "SNX9", "TNFRSF18", "ZBED2", "TNFRSF9", "VCAM1",
                         "CXCL13", "GEM")
    yost_module <- yost_exhaustion[yost_exhaustion %in% var_feats]
    integrated_gex <- AddModuleScore(integrated_gex,
                                     features = list(yost_module),
                                     name = "Yost_exhaustion_score")
    # markers from Yost et al.
    yost_markers <- c("CD3D","CD8A", "CD4", "FOXP3", "CCR7", "IL26", "CD200",
                     "EOMES", "KLRD1", "IFNG", "HAVCR2", "CXCL13", "ENTPD1", "PDCD1", "ITGAE")
    yost_markers <- yost_markers[yost_markers %in% var_feats]
    yost_fp <- FeaturePlot(integrated_gex,
                           features = c(yost_markers,
                                        "Yost_exhaustion_score1"))
    ggsave(filename = paste0(in_list$save_dir, "/GEX_UMAP_YostFP.png"),
           plot = yost_fp,
           width = 12,
           height = 8)
    # set of genes identified from literature to be plotted individually
    plot_dir <- paste0(in_list$save_dir, "/Feature_plots_GEX/")
    if (!dir.exists(plot_dir)) {
      dir.create(plot_dir)
    }
    suppressWarnings(gen_feature_plots(gex_seurat = integrated_gex,
                                       plot_dir = plot_dir))
  }

  return(NULL)
}
