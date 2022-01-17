VDJTOOLS="java -Xmx20G -jar /cygene/software/biosoftware/vdjtools/vdjtools-1.2.1/vdjtools-1.2.1.jar"

### basic statistics:
## calc Basic STats , read counts, avg clonotype, nonfunctional clonotype,
VDJTOOLS CalcBasicStats -m data/metadata.txt out/0

###
VDJTOOLS CalcSegmentUsage -m  data/metadata.txt -p -f age -n out/2

# calculate spectratype : CDR3 read histogram:
VDJTOOLS CalcSpectratype -m data/metadata.txt out/1	

# VDJTOOLS
VDJTOOLS PlotFancySpectratype data/A4-i125.txt.gz out/3

# 
VDJTOOLS PlotFancyVJUsage data/A4-i125.txt.gz out/5

# 
VDJTOOLS PlotSpectratypeV data/A4-i125.txt.gz out/4

# diversity
VDJTOOLS PlotQuantileStats data/A4-i125.txt.gz out/6

# Rarefaction Plot
VDJTOOLS RarefactionPlot -m data/metadata.txt -f age -n -l sample.id out/8

VDJTOOLS CalcDiversityStats -m data/metadata.txt out/7

# overlap between two samples
VDJTOOLS OverlapPair -p data/A4-i189.txt.gz data/A4-i190.txt.gz out/9

# 
VDJTOOLS CalcPairwiseDistances -m data/metadata.small.txt out/10

#
VDJTOOLS TrackClonotypes -m data/metadata.small.txt -f age -x 0 -p out/11

java -Xmx20G -jar /cygene/software/biosoftware/vdjtools/vdjtools-1.2.1/vdjtools-1.2.1.jar CalcBasicStats /cygene2/work/P0000-Blackbird/2103-BL001/BL001004/BL14-PBMCs-CD3-PD1orCD39.T.cells_G328E2/G328E2L1_TCRseq_IMGT/outs/clonotypes.csv out/h
java -Xmx20G -jar /cygene/software/biosoftware/vdjtools/vdjtools-1.2.1/vdjtools-1.2.1.jar Convert /cygene2/work/P0000-Blackbird/2103-BL001/BL001004/BL14-PBMCs-CD3-PD1orCD39.T.cells_G328E2/G328E2L1_TCRseq_IMGT/outs/clonotypes.csv out/convert

