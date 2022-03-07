# Enter commands in R (or R studio, if installed)
install.packages('Seurat')
library(Seurat)

packageVersion("Seurat")

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("GenVisR")
library(GenVisR)

BiocManager::install("TCGAbiolinks")
library(TCGAbiolinks)

?GDCquery_Maf
maf_file = GDCquery_Maf("COAD",save.csv = F, directory = "GDCdata", pipelines = "mutect2")
maf_file = GDCquery_Maf("COAD",save.csv = F, directory = "GDCdata", pipelines = "muse") # muse, somaticsniper,varscan2

maf_file[1:15,c("Tumor_Sample_Barcode", "Hugo_Symbol", "Variant_Classification")]

my_maf = read.csv("/cygene2/work/P0000-Blackbird/2103-CR001/CR001003/CR001003_WES_RNA/CR001003.variants.funcotated.without.header.MAF.xls",header = T, sep = "\t")

my_maf[1:15,c("Tumor_Sample_Barcode", "Hugo_Symbol", "Variant_Classification")]

waterfall(maf_file, fileType = "MAF")
factor(maf_file[,"Variant_Classification"])
length(which(maf_file$Variant_Classification =="Splice_Region"))
maf_file = subset(maf_file, Variant_Classification != "Splice_Region")
waterfall(maf_file, fileType = "MAF")
waterfall(my_maf, fileType = "MAF", mainRecurCutoff = 0.06, top=20)
########################### 2
# Create input data
data <- brcaMAF[brcaMAF$Hugo_Symbol == "TP53", c("Hugo_Symbol", "amino_acid_change_WU")]
head(data)
data <- as.data.frame(cbind(data, "ENST00000269305"))
colnames(data) <- c("gene", "amino_acid_change", "transcript_name")
data
# Call lolliplot
lolliplot(data)
#### my data
maf_file[1:10,]
library(dplyr)
data = maf_file %>% select(Hugo_Symbol, HGVSp_Short, Transcript_ID) %>%
  subset(Hugo_Symbol == "TP53") %>% filter(!is.na(HGVSp_Short)) %>%
  rename(gene = Hugo_Symbol, amino_acid_change = HGVSp_Short, transcript_name=Transcript_ID)
data = my_maf %>% select(Hugo_Symbol, Protein_Change, Annotation_Transcript) %>% 
  subset(Hugo_Symbol =="TP53") %>% filter(!is.na(Protein_Change)) %>%
  rename(gene = Hugo_Symbol, amino_acid_change = Protein_Change, transcript_name=Annotation_Transcript) 

head(data)
library(stringr)
data$transcript_name = str_replace_all(data$transcript_name,".8","")
lolliplot(data)
##################################### 3 
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(BSgenome.Hsapiens.UCSC.hg38)
library(BSgenome.Hsapiens.UCSC.hg19)
BiocManager::install(c("TxDb.Hsapiens.UCSC.hg38.knownGene","TxDb.Hsapiens.UCSC.hg19.knownGene","BSgenome.Hsapiens.UCSC.hg19","BSgenome.Hsapiens.UCSC.hg38"))

txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
genome <- BSgenome.Hsapiens.UCSC.hg19

gr <- GRanges(seqnames = c("chr10"), ranges=IRanges(start = c(89622195), end = c(89729532)),strand = strand(c("+")))
# Create Data for input
start <- c(89622194:89729524)
end <- c(89622195:89729525)
chr <- 10
cov <- c(rnorm(1e+05, mean = 40), rnorm(7331, mean = 10))
cov_input_A <- as.data.frame(cbind(chr, start, end, cov))

start <- c(89622194:89729524)
end <- c(89622195:89729525)
chr <- 10
cov <- c(rnorm(50000, mean = 40), rnorm(7331, mean = 10), rnorm(50000, mean = 40))
cov_input_B <- as.data.frame(cbind(chr, start, end, cov))

data <- list(`Sample A` = cov_input_A, `Sample B` = cov_input_B)
data
# Call genCov
genCov(data, txdb, gr, genome, gene_labelTranscriptSize = 2, transform = NULL, base = NULL)
genCov(data, txdb, gr, genome, transform = c("Intron", "CDS", "UTR"), base = c(10, 2, 2), reduce = TRUE)
################################### 4 TvTi
p1 = TvTi(maf_file, lab_txtAngle = 90, fileType = "MAF")
p2 = TvTi(my_maf, fileType = "MAF")
p1 | p2

#################################### 5 cnSpec
LucCNseg
cnSpec(LucCNseg, genome = "hg38")
cnSpec(LucCNseg, y=hg19chr)

##################################### 
# Create data
chromosome <- "chr14"
coordinate <- sort(sample(0:106455000, size = 2000, replace = FALSE))
cn <- c(rnorm(300, mean = 3, sd = 0.2), rnorm(700, mean = 2, sd = 0.2), rnorm(1000, mean = 3, sd = 0.2))
data <- as.data.frame(cbind(chromosome, coordinate, cn))
data
# Call cnView with basic input
cnView(data,              chr = "chr14", genome = "hg19", ideogram_txtSize = 4)

# create segment data
dataSeg <- data.frame(chromosome = c(14, 14, 14), start = coordinate[c(1, 301, 1001)],
                      end = coordinate[c(300, 1000, 2000)], segmean = c(3, 2, 3))
dataSeg
cnView(data, z = dataSeg, chr = "chr14", genome = "hg19", ideogram_txtSize = 4)

######################################## 7
# Example input to x
x <- matrix(sample(1e+05, 500), nrow = 50, ncol = 10, dimnames = list(0:49, paste0("Sample",1:10)))

covBars(x)
######################################## 8 
cnFreq(LucCNseg)
head(LucCNseg)
######################################## 9
# Obtain cytogenetic information for the genome of interest
data <- cytoGeno[cytoGeno$genome == "hg38", ]
head(data)
# Call ideoView for chromosome 1
ideoView(data, chromosome = "chr1", txtSize = 4)
######################################## 10
head(HCC1395_Germline)
lohSpec(x = HCC1395_Germline)
