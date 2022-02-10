# Basic analysis
# java -jar ../vdjtools-1.2.1.jar Convert -S mixcr sample1.mixcr.TRA.txt vdjtools_Format
java -jar ../vdjtools-1.2.1.jar Convert -S mixcr -m metadata.txt  vdjtools_Format


##################
##################
## Basic analysis
# CalcBasicStats
java -jar ../vdjtools-1.2.1.jar CalcBasicStats   -m metadata.txt basic_01_CalcBasicStats
java -jar ../vdjtools-1.2.1.jar CalcSegmentUsage -m metadata.txt --plot --plot-type png basic_02_CalcSegmentUsage
java -jar ../vdjtools-1.2.1.jar CalcSpectratype  -m metadata.txt basic_03_CalcSpectratype_nt
java -jar ../vdjtools-1.2.1.jar CalcSpectratype  -m metadata.txt --amino-acid basic_03_CalcSpectratype_aa

# plot a histogram pdf 
java -jar ../vdjtools-1.2.1.jar PlotFancySpectratype  --top 10 --plot-type png vdjtools_Format.sample1.txt basic_04_PlotFancySpectratype_sample1
java -jar ../vdjtools-1.2.1.jar PlotFancySpectratype  --top 10 --plot-type png vdjtools_Format.sample2.txt basic_04_PlotFancySpectratype_sample2
java -jar ../vdjtools-1.2.1.jar PlotFancySpectratype  --top 10 --plot-type png vdjtools_Format.sample3.txt basic_04_PlotFancySpectratype_sample3

# plot a circos-style V-J usage plot Error with Rscript, 
java -jar ../vdjtools-1.2.1.jar PlotFancyVJUsage --plot-type png vdjtools_Format.sample1.txt basic_05_circos_sample1
java -jar ../vdjtools-1.2.1.jar PlotFancyVJUsage --plot-type png vdjtools_Format.sample2.txt basic_05_circos_sample2
java -jar ../vdjtools-1.2.1.jar PlotFancyVJUsage --plot-type png vdjtools_Format.sample3.txt basic_05_circos_sample3

# plot a detailed spectratype containing additional info displays CDR3 length distribution for clonotypes from top N Variable segment families. This plot is useful to detect type 1 and type 2 repertoire biases, that could arise under pathological conditions.
java -jar ../vdjtools-1.2.1.jar PlotSpectratypeV --top 12 --plot-type png vdjtools_Format.sample1.txt basic_06_PlotSpectratypeV_sample1
java -jar ../vdjtools-1.2.1.jar PlotSpectratypeV --top 12 --plot-type png vdjtools_Format.sample2.txt basic_06_PlotSpectratypeV_sample2
java -jar ../vdjtools-1.2.1.jar PlotSpectratypeV --top 12 --plot-type png vdjtools_Format.sample3.txt basic_06_PlotSpectratypeV_sample3


###################
###################  Part 2.
# Diversity estimation
# Plots a three-layer donut chart to visualize the repertoire clonality.
java -jar ../vdjtools-1.2.1.jar PlotQuantileStats --top 5 --plot-type png vdjtools_Format.sample1.txt Diversity_est_01_donut_chart_sample1
java -jar ../vdjtools-1.2.1.jar PlotQuantileStats --top 5 --plot-type png vdjtools_Format.sample2.txt Diversity_est_01_donut_chart_sample2
java -jar ../vdjtools-1.2.1.jar PlotQuantileStats --top 5 --plot-type png vdjtools_Format.sample3.txt Diversity_est_01_donut_chart_sample3

# Plots rarefaction curves for specified list of samples, that is, the dependencies between sample diversity and sample size.
java -jar ../vdjtools-1.2.1.jar RarefactionPlot -m metadata.txt --plot-type png Diversity_est_02_rarefaction_curve

# Computes a set of diversity statistics,  None Graphical output
java -jar ../vdjtools-1.2.1.jar CalcDiversityStats -m metadata.txt  Diversity_est_03_DiversityStats


##################
##################
# Repertoire overlap analysis
# overlapPair
java -jar ../vdjtools-1.2.1.jar OverlapPair --top 20 --plot --plot-type png vdjtools_Format.sample1.txt vdjtools_Format.sample2.txt  Overlap_01_overlapPair_sample1_sample2
java -jar ../vdjtools-1.2.1.jar CalcPairwiseDistances  --plot --plot-type png -m metadata.txt CalcPairwiseDistances_02

# multiple sample
java -jar ../vdjtools-1.2.1.jar ClusterSamples --plot --plot-type png CalcPairwiseDistances_02 ClusterSamples_03
java -jar ../vdjtools-1.2.1.jar TestClusters   --plot-type png ClusterSamples_03 TestClusters_04

java -jar ../vdjtools-1.2.1.jar TrackClonotypes -p --plot-type png -m metadata.txt --top 100 TrackClonotypes_05

#################
#################
# Correct
java -jar ../vdjtools-1.2.1.jar  G19E6L2.mixcr.out.clonotypes.TRB.txt 041
# Decontaminate
java -jar ../vdjtools-1.2.1.jar  G19E6L2.mixcr.out.clonotypes.TRB.txt 041
# DownSample
java -jar ../vdjtools-1.2.1.jar DownSample  G19E6L2.mixcr.out.clonotypes.TRB.txt 041
# FilterNonFunctional
java -jar ../vdjtools-1.2.1.jar FilterNonFunctional  G19E6L2.mixcr.out.clonotypes.TRB.txt 041
# SelectTop
java -jar ../vdjtools-1.2.1.jar SelectTop  G19E6L2.mixcr.out.clonotypes.TRB.txt 041
# FilterByFrequency
java -jar ../vdjtools-1.2.1.jar FilterByFrequency  G19E6L2.mixcr.out.clonotypes.TRB.txt 041
# ApplySampleAsFilter
java -jar ../vdjtools-1.2.1.jar ApplySampleAsFilter  G19E6L2.mixcr.out.clonotypes.TRB.txt 041
# FilterBySegment
java -jar ../vdjtools-1.2.1.jar FilterBySegment  G19E6L2.mixcr.out.clonotypes.TRB.txt 041

################
################
# Operate on clonotype table
# JoinSamples
java -jar ../vdjtools-1.2.1.jar JoinSamples --plot G19E6L2.mixcr.out.clonotypes.TRB.txt 051
# PoolSamples
java -jar ../vdjtools-1.2.1.jar PoolSamples --plot G19E6L2.mixcr.out.clonotypes.TRB.txt 051

###############
###############
## Annotation
# CalcDegreeStats
java -jar ../vdjtools-1.2.1.jar CalcDegreeStats  G19E6L2.mixcr.out.clonotypes.TRB.txt 061
# CalcCdrAAProfile
java -jar ../vdjtools-1.2.1.jar CalcCdrAAProfile  G19E6L2.mixcr.out.clonotypes.TRB.txt 062
# Annotate
java -jar ../vdjtools-1.2.1.jar Annotate  G19E6L2.mixcr.out.clonotypes.TRB.txt 063
# ScanDatabase(VDJmatch)
java -jar /cygene/software/biosoftware/vdjmatch-1.3.1/vdjmatch-1.3.1.jar  match G19E6L2.mixcr.out.clonotypes.TRB.txt 064
java -jar /cygene/software/biosoftware/vdjmatch-1.3.1/vdjmatch-1.3.1.jar  cluster G19E6L2.mixcr.out.clonotypes.TRB.txt 065


###############
###############
## Utilities
# SplitMetadata
java -jar ../vdjtools-1.2.1.jar SplitMetadata metadata.txt output_dir
# FilterMetadata
java -jar ../vdjtools-1.2.1.jar FilterMetadata metadata.txt output_dir output_suffix 
# Convert
java -jar ../vdjtools-1.2.1.jar Convert -S mixcr G19E6L2.mixcr.out.clonotypes.TRB.txt output_prefix
# RInstall
java -jar ../vdjtools-1.2.1.jar RInstall


