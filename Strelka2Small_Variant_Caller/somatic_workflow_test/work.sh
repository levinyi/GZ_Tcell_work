/cygene/software/biosoftware/manta/bin/configManta.py \
--normalBam /cygene2/work/P0000-Blackbird/2103-BL001/BL001004/BL001004_WES_RNAseq/BL001004/BL001004.Normal.duplicates_marked_sorted_fixed.BQSR.bam \
--tumorBam /cygene2/work/P0000-Blackbird/2103-BL001/BL001004/BL001004_WES_RNAseq/BL001004/BL001004.Tumor.duplicates_marked_sorted_fixed.BQSR.bam \
--referenceFasta /cygene/database/GATK_resource_bundle/hg38/Homo_sapiens_assembly38.fasta \
--runDir ./manta_results

./manta_results/runWorkflow.py --quiet


# configuration
STRELKA_INSTALL_PATH="/cygene/software/biosoftware/strelka/strelka-2.9.10.centos6_x86_64/"
MANTA_ANALYSIS_PATH="./manta_results"
${STRELKA_INSTALL_PATH}/bin/configureStrelkaSomaticWorkflow.py \
    --normalBam /cygene2/work/P0000-Blackbird/2103-BL001/BL001004/BL001004_WES_RNAseq/BL001004/BL001004.Normal.duplicates_marked_sorted_fixed.BQSR.bam \
    --tumorBam /cygene2/work/P0000-Blackbird/2103-BL001/BL001004/BL001004_WES_RNAseq/BL001004/BL001004.Tumor.duplicates_marked_sorted_fixed.BQSR.bam \
    --referenceFasta /cygene/database/GATK_resource_bundle/hg38/Homo_sapiens_assembly38.fasta \
    --runDir strelka_somatic_result \
    --indelCandidates ${MANTA_ANALYSIS_PATH}/results/variants/candidateSmallIndels.vcf.gz \
    --exome \
    --callRegions Regions_hg38.bed.gz 
# execution on a single local machine with 72 parallel jobs
demo_somatic/runWorkflow.py --quiet -m local -j 72
