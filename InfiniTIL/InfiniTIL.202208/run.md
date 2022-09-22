# $10X is the absolute path to the base dir of this pipeline and data
# modify 01.TCR_selection.R:
# line 7:  set the correct path to InfiniTIL
# line 50: set the correct path to InfiniTIL/index/science.abl5447_table_s4.xlsx

# step 1: post-processing of CellRanger output
cd $10X/InfiniTIL
Rscript BatchSeurat.R -i $10X/data/HC29
# manual QC
# some commonly used thresholds:
# pct_mito < 10%; pct_ribo > 10%;

# step 2: calculate AUC scores for NeoTCR4 and NeoTCR8
cd $10X/data/HC29/analysis
Rscript /cygene2/pipeline/10X/InfiniTIL.202208/InfiniTIL.202208/01.TCR_selection.R

# step 3: manually determine cluster identities (CD4, CD4_Treg, and CD8)
# see GEX_UMAP_GexClusterByRes.png for clustering results
# see GEX_UMAP_YostFP.png for expression patterns of CD4 (for CD4), FOXP3 (for CD4_Tref), and CD8A (for CD8) 

# step 3: prepare data for InfiniTIL analysis
# remember to modify cluster membership information on lines 40-42 in $10X/InfiniTIL/02.get_InfiniTIL_data.R
ln -s /cygene2/pipeline/10X/InfiniTIL.202208/10X-analysis.202208/src .
Rscript /cygene2/pipeline/10X/InfiniTIL.202208/InfiniTIL.202208/02.get_InfiniTIL_data.R

# step 4: run InfiniTIL
perl /cygene2/pipeline/10X/InfiniTIL.202208/InfiniTIL.202208/03.InfiniTIL.pl tcr_stats_barcodes.csv merged_tcr_table.csv InfiniTIL_data.csv

# step 5: manual selection of CB

# step 6: prepare data for CD4 and CD8 plots
perl /cygene2/pipeline/10X/InfiniTIL.202208/InfiniTIL.202208/04.selected_CB.pl tcr_stats_barcodes.csv InfiniTIL_data.csv selected_convCD4.csv TCR4 > convCD4.out
perl /cygene2/pipeline/10X/InfiniTIL.202208/InfiniTIL.202208/04.selected_CB.pl tcr_stats_barcodes.csv InfiniTIL_data.csv selected_TregCD4.csv TCR4 > TregCD4.out
perl /cygene2/pipeline/10X/InfiniTIL.202208/InfiniTIL.202208/04.selected_CB.pl tcr_stats_barcodes.csv InfiniTIL_data.csv selected_CD8.csv TCR8 > CD8.out

# step 7: make the plot
Rscript /cygene2/pipeline/10X/InfiniTIL.202208/InfiniTIL.202208/05.InfiniTIL_plot.R
