# File "barcodes.tsv" comes from: /cygene2/work/P0000-Blackbird/2103-GB001/GB001003/analysis/GB001003_G172E2L1_RNAseq/outs/filtered_feature_bc_matrix
/cygene/software/biosoftware/scHLAcount-0.1.0/sc_hla_count_linux \
	--bam /cygene2/work/P0000-Blackbird/2103-GB001/GB001003/analysis/GB001003_G172E2L1_RNAseq/outs/possorted_genome_bam.bam  \
	--cell-barcodes barcodes.tsv \
	-d /cygene/database/imgt-hla

