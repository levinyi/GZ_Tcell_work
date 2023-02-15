#' scarHRD is an R package which determines the levels of homologous recombination deficiency (telomeric allelic imbalance, 
#' loss off heterozygosity, number of large-scale transitions) based on NGS (WES, WGS) data.
#' The first genomic scar based homologous recombination deficiency measures were produced using SNP arrays. 
#' Since this technology has been largely replaced by next generation sequencing it has become important to develop algorithms 
#' that derive the same type of genomic scar-scores from next generation sequencing (WXS, WGS) data. In order to perform this 
#' analysis, here we introduce the scarHRD R package and show that using this method the SNP-array based and next generation 
#' sequencing based derivation of HRD scores show good correlation.


# A modification of the copynumber R package needs to be used which can be installed via devtools from github:
library(devtools)
install_github('aroneklund/copynumber')

# install_bitbucket('sequenza_tools/sequenza') 
install.packages("sequenza")


install_github('sztup/scarHRD',build_vignettes = TRUE)


library("scarHRD")
scar_score("/cygene2/work/IIT-WES-RNAseq/SA501001/SA501001_NeoAg_ID/Sequenza-analysis/SA501001_TU-gDNA.small.seqz.gz", reference = "grch38", seqz=TRUE)

scar_score("F:/Documents/scarHRD/examples/test2.txt",reference = "grch38", seqz=FALSE)

#' \seg -- input file name \reference -- the reference genome used, \grch38 or \grch37 or \mouse (default: grch38) 
#' \seqz -- \TRUE if the input file is a smallo.seqz.gz file, otherwise \FALSE (default: TRUE) 
#' \ploidy -- optional, previously estimated ploidy of the sample 
#' \outputdir -- optional, the path to the output directory 
#' \chr.in.names -- optional, default: TRUE, set to FALSE if input file does not contain 'chr' in chromosome names.

ref: https://rdrr.io/github/sztup/scarHRD/f/README.md
