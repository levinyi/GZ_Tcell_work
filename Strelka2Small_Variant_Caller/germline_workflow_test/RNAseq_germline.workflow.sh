# configuration
STRELKA_INSTALL_PATH="/cygene/software/biosoftware/strelka/strelka-2.9.10.centos6_x86_64/"
MANTA_ANALYSIS_PATH="."
${STRELKA_INSTALL_PATH}/bin/configureStrelkaGermlineWorkflow.py \
    --bam /cygene2/work/P0000-Blackbird/2103-BL001/BL001004/BL001004_WES_RNAseq/BL001004/BL001004.RNA.Aligned.sortedByCoord.out.bam \
    --referenceFasta /cygene/database/GATK_resource_bundle/hg38/Homo_sapiens_assembly38.fasta \
    --runDir RNA_germline

# execution on a single local machine with 72 parallel jobs
./RNA_germline/runWorkflow.py -m local -j 72 --quiet

