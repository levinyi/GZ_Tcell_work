
### -----------------------
### This is a Shiny app for performing QC on sets of 10X data using Seurat
### Author: Seon Kinrot (RootPath Genomics)
### -----------------------
### requirements
### -----------------------
library(shiny)
library(dplyr)
library(ggplot2)
library(patchwork)
library(here)
source(here("src", "SeuratWrappers.R"))
### -----------------------
### load data and set up global variables
### -----------------------
data_dir <- here("ShinyQC", "data")
input_fl <- paste0(data_dir, "/input.rds")
output_fl <- paste0(data_dir, "/output.rds")
if (!file.exists(input_fl)) {
  cat("No input file detected!\n")
  return(NULL)
}
# data should be a list of Seurat objects
data <- readRDS(input_fl)
projects <- names(data)
if (is.null(projects)) {
  projects <- as.character(seq_along(data))
  names(data) <- projects
}
# initialize some output and default filters
data_out <- data
default_filters <- list("mito_th"=c(0, 100),
                        "ribo_th"=c(0, 100),
                        "gene_th"=c(0, Inf),
                        "umi_th"=c(0, Inf))
filter_list <- vector(mode="list", length=length(projects))
names(filter_list) <- projects
for (project in projects) {
  filter_list[[project]] <- default_filters
}

### -----------------------
### UI
### -----------------------
ui <- fluidPage(

  ## title
  titlePanel(textOutput("title")),

  ## page layout is side + main panels
  sidebarLayout(
    ## sidebar
    sidebarPanel(

      # species selection radio panel
      radioButtons("species",
                   label="Species",
                   choices=c("human", "mouse")),
      br(),

      # dataset selection dropdown tool
      selectInput("dataset",
                  label="Dataset to QC",
                  choices=projects),

      # horizontal ruler to separate inputs from summaries
      hr(),

      # help text with summary statistics
      textOutput("text_summary"),
      br(),

      # download/ save button for data_out
      # note that the program also auto-saves upon session end
      downloadButton("save_data", "Download post-QC data")
    ),  # end sidebarPanel

    ## main panel
    mainPanel(

      # row with all the interactive UIs
      fluidRow(
        column(4, uiOutput("gene_ui")),
        column(4, uiOutput("mito_ui")),
        column(4, uiOutput("ribo_ui"))
      ),

      hr(),

      # buttons for updating plots, resetting etc.
      fluidRow(
        column(4, actionButton("reset", "Reset filters")),
        column(4, actionButton("apply_filter", "Apply filters")),
        column(4, downloadButton("save_plot", "Download QC plots"))
      ),
      br(), br()
    )  # end mainPanel
  )  # end sidebarLayout
)  # end fluidPage/ ui

### -----------------------
### server
### -----------------------
server <- function(input, output, session) {
  ## title text
  output$title <- renderText({
    paste("Performing QC for", input$dataset)
  })

  ## helper functions to save the output data
  # accepts a file to comply with downloadHandler format, but can be called
  # with no arguments to save to output_fl
  download_data <- function(file=NULL) {
    if (is.null(file)) {
      file <- output_fl
    }
    saveRDS(data_out, file)
  }

  ## download handler for data download
  output$save_data <- downloadHandler(
    filename = function() {
      output_fl
    },
    content = download_data
  )

  ## reactive objects for:
  # selected species
  sample_species <- reactive({
    input$species
  })
  # name of current dataset
  sample_name <- reactive({
    input$dataset
  })
  # reactive filters from input (may not be applied)
  cur_filters <- reactive({
    list("mito_th"=input$mito_th,
         "ribo_th"=input$ribo_th,
         "gene_th"=input$gene_th,
         "umi_th"=input$umi_th)
  })
  # max allowed values for gene and umi count filters
  max_genes <- reactive({
    gene_counts <- FetchData(object = data[[input$dataset]],
                             vars = "nFeature_RNA")
    ceiling(log10(max(gene_counts, na.rm=TRUE)))
  })
  max_umis <- reactive({
    umi_counts <- FetchData(object = data[[input$dataset]],
                            vars = "nCount_RNA")
    ceiling(log10(max(umi_counts, na.rm=TRUE)))
  })

  ## helper function - updates saved filters to curernt slider input values
  apply_filters <- function() {
    filter_list[[sample_name()]] <<- cur_filters()
  }

  ## helper function that just runs QC within a tryCatch
  shiny_qc <- function(object, filters, species) {
    qc_out <- suppressWarnings(
      tryCatch(
        {run_qc(srt_obj = object,
                min_fsets = c(filters$mito_th[1], filters$ribo_th[1]),
                max_fsets = c(filters$mito_th[2], filters$ribo_th[2]),
                min_feats = 10**(filters$gene_th[1]),
                max_feats = 10**(filters$gene_th[2]),
                min_umi = 10**(filters$umi_th[1]),
                max_umi = 10**(filters$umi_th[2]),
                species = species,
                return_plot = TRUE, return_object = TRUE)
        },  # end main clause, start error clause
        error = function(cond) {
          message(paste0("Failed to run normally for sample ", sample_name(),
                         "; skipping CC scoring"))
          run_qc(srt_obj = object,
                 min_fsets = c(filters$mito_th[1], filters$ribo_th[1]),
                 max_fsets = c(filters$mito_th[2], filters$ribo_th[2]),
                 min_feats = 10**(filters$gene_th[1]),
                 max_feats = 10**(filters$gene_th[2]),
                 min_umi = 10**(filters$umi_th[1]),
                 max_umi = 10**(filters$umi_th[2]),
                 species = species, score_cc = FALSE,
                 return_plot = TRUE, return_object = TRUE)
        }  # end error clause
      )  # end tryCatch
    )  # end suppressWarnings
    return(qc_out)
  }  # end shiny_qc

  ## helper function - renders the current plots
  render_plots <- function(qc_plots) {
    flt <- cur_filters()
    output$umi_gene_scatter <- renderPlot({
      qc_plots[["umi_gene_scatter"]]
    })
    output$umi_mito_scatter <- renderPlot({
      qc_plots[["umi_pct_mito_scatter"]]
    })
    output$umi_ribo_scatter <- renderPlot({
      qc_plots[["umi_pct_ribo_scatter"]]
    })
    output$gene_hist <- renderPlot({
      qc_plots[["gene_hist"]]
    })
    output$mito_hist <- renderPlot({
      qc_plots[["pct_mito_hist"]]
    })
    output$ribo_hist <- renderPlot({
      qc_plots[["pct_ribo_hist"]]
    })
  }

  ## helper function - update plots and the filtered current dataset
  # returns a list of plots
  update_plots <- function(object, filters, species) {
    # run QC
    qc_out <- shiny_qc(object, filters, species)
    # render plots
    render_plots(qc_out$plots)
    # render text summary
    n_cells <- c(dim(qc_out$filtered_seurat)[2],  # number of barcodes post-QC
                 dim(object)[2])  # number of barcodes pre-QC
    output$text_summary <- renderText({
      paste0("Cells retained after filter:\n",
             n_cells[1], " out of ", n_cells[2],
             " (", round(100*n_cells[1]/n_cells[2]), "%)")
    })  # end renderText
    return(qc_out)
  }  # end update_plots

  ## helper function - saves filtered dataset to data_out
  save_srt <- function(filtered_obj) {
    data_out[[sample_name()]] <<- filtered_obj
  }

  ## when user clicks on "apply filters", save slider inputs
  ## then update plots and data_out
  observe({
    apply_filters()
    qc_out <- update_plots(object = data[[sample_name()]],
                           filters = filter_list[[sample_name()]],
                           species = sample_species())
    save_srt(qc_out$filtered_seurat)
  }) %>%
    bindEvent(input$apply_filter)

  ## helper function - updates slider inputs while considering limitations
  update_sliders <- function(default=FALSE) {
    if (default) {
      flt <- default_filters
    }
    else {
      flt <- filter_list[[sample_name()]]
    }
    # enforce max allowed values for:
    # genes
    mg <- max_genes()
    if (flt$gene_th[2] > mg) {
      flt$gene_th[2] <- mg
    }
    # UMIs
    mu <- max_umis()
    if (flt$umi_th[2] > mu) {
      flt$umi_th[2] <- mu
    }
    # update sliders
    updateSliderInput(session, "gene_th", value = flt$gene_th)
    updateSliderInput(session, "umi_th", value = flt$umi_th)
    updateSliderInput(session, "mito_th", value = flt$mito_th)
    updateSliderInput(session, "ribo_th", value = flt$ribo_th)
  }

  ## when changing datasets, load saved filters and update sliders
  ## then redo plots (loads correct Seurat object internally)
  observe({
    update_sliders()
    update_plots(object = data[[sample_name()]],
                 filters = filter_list[[sample_name()]],
                 species = sample_species())
  }) %>%
    bindEvent(input$dataset)

  ## when reset button is clicked, update sliders to their default values
  ## then apply these filters
  observe({
    update_sliders(default=TRUE)
    apply_filters()
    qc_out <- update_plots(object = data[[sample_name()]],
                           filters = filter_list[[sample_name()]],
                           species = sample_species())
    save_srt(qc_out$filtered_seurat)
  }) %>%
    bindEvent(input$reset)

  ## helper function to save plots
  download_plot <- function(file) {
    qc_out <- update_plots(object = data[[sample_name()]],
                           filters = filter_list[[sample_name()]],
                           species = sample_species())
    qc_plots <- qc_out$plots
    combined_plot <- qc_plots[[1]]
    for (iplot in seq(2, length(qc_plots))) {
      combined_plot <- combined_plot + qc_plots[[iplot]]
    }
    ggsave(file, plot = combined_plot, width = 8, height = 6)
  }

  ## download handler for plot download
  output$save_plot <- downloadHandler(
    filename = function() {
      paste0(data_dir, "/", input$dataset, "_QCplots.png")
    },
    content = download_plot
  )

  ## reactive plot UIs for each plot:
  # plots of gene counts
  output$gene_ui <- renderUI({
    # initiate filters
    flt <- filter_list[[input$dataset]]
    vals_g <- flt$gene_th
    vals_u <- flt$umi_th
    # enforce max allowed thresholds for:
    # genes
    mg <- max_genes()
    if (vals_g[2] > mg) {
      vals_g[2] <- mg
    }
    # UMIs
    mu <- max_umis()
    if (vals_u[2] > mu) {
      vals_u[2] <- mu
    }

    # column object to be returned
    column(12,
          # scatter plot
          fluidRow(plotOutput("umi_gene_scatter")),
          hr(),
          # sliders for filter thresholds
          fluidRow(
            wellPanel(
              sliderInput("gene_th", "Threshold (log10 gene count):",
                          min = 0, max = mg, step = 0.05, value = vals_g),
              sliderInput("umi_th", "Threshold (log10 UMI count):",
                          min = 0, max = mu, step = 0.05, value = vals_u)
            )  # end wellPanel
          ),  # end fluidRow
          hr(),
          # histogram
          fluidRow(plotOutput("gene_hist"))
    )  # end column
  })  # end gene_ui

  # plots of Mitochondrial gene content
  output$mito_ui <- renderUI({
    column(12,
          # scatter plot
          fluidRow(plotOutput("umi_mito_scatter")),
          hr(),
          # sliders for filter thresholds
          fluidRow(
            wellPanel(
              sliderInput("mito_th", "Threshold (%Mitochondrial):",
                          min = 0, max = 100, step = 0.5,
                          value = filter_list[[input$dataset]]$mito_th),
              br(), br(), br(), br(), br(), br()
            )  # end wellPanel
          ),  # end fluidRow
          hr(),
          # histogram
          fluidRow(plotOutput("mito_hist"))
    )  # end column
  })  # end mito_ui

  # plots for Ribosomal gene content
  output$ribo_ui <- renderUI({
    column(12,
          # scatter plot
          fluidRow(plotOutput("umi_ribo_scatter")),
          hr(),
          # sliders for filter thresholds
          fluidRow(
            wellPanel(
              sliderInput("ribo_th", "Threshold (%Ribosomal):",
                          min = 0, max = 100, step = 0.5,
                          value = filter_list[[input$dataset]]$ribo_th),
              br(), br(), br(), br(), br(), br()
            )  # end wellPanel
          ),  # end fluidRow
          hr(),
          # histogram
          fluidRow(plotOutput("ribo_hist"))
    )  # end column
  })  # end ribo_ui

  ## actions taken when session ends
  session$onSessionEnded(function() {
    cat("Saving post-QC data to disk...\n")
    download_data()
    stopApp()
  })
}  # end server

### -----------------------
### app call
shinyApp(ui, server)
