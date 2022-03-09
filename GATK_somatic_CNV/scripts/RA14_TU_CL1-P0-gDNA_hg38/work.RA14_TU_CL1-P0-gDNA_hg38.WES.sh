/cygene/software/biosoftware/bwa-0.7.17/bwa mem -t 30 -Y -H "@HD\tVN:1.5\tGO:none\tSO:coordinate" -R "@RG\tID:RA14_TU_CL1-P0-gDNA_Normal\tSM:Normal\tLB:WES\tPL:ILLumina\tPU:HVW2MCCXX:6:none" /cygene/database/GATK_resource_bundle/hg38/Homo_sapiens_assembly38.fasta /cygene2/data/2020/Blackbird/2nd-batch-5samples/20200624_H5JGWDSXY-X101SC20031639-Z01_Result/Rawdata/RA001004-BC-gDNA_L3_1.fq.gz /cygene2/data/2020/Blackbird/2nd-batch-5samples/20200624_H5JGWDSXY-X101SC20031639-Z01_Result/Rawdata/RA001004-BC-gDNA_L3_2.fq.gz | /usr/bin/samtools view -@ 30 -buhS -t /cygene/database/GATK_resource_bundle/hg38/Homo_sapiens_assembly38.fasta.fai - | /usr/bin/samtools sort -@ 30 -o RA14_TU_CL1-P0-gDNA.Normal.sortedByCoord.bam - 
/cygene/software/biosoftware/bwa-0.7.17/bwa mem -t 30 -Y -H "@HD\tVN:1.5\tGO:none\tSO:coordinate" -R "@RG\tID:RA14_TU_CL1-P0-gDNA_Tumor\tSM:Tumor\tLB:WES\tPL:ILLumina\tPU:HVW2MCCXX:6:none" /cygene/database/GATK_resource_bundle/hg38/Homo_sapiens_assembly38.fasta /cygene2/data/2022/20220207_SF-X101SC20031639-Z01_Result/Rawdata/RA14_TU_CL1-P0-gDNA_2022112_L1_1.fq.gz /cygene2/data/2022/20220207_SF-X101SC20031639-Z01_Result/Rawdata/RA14_TU_CL1-P0-gDNA_2022112_L1_2.fq.gz | /usr/bin/samtools view -@ 30 -buhS -t /cygene/database/GATK_resource_bundle/hg38/Homo_sapiens_assembly38.fasta.fai - | /usr/bin/samtools sort -@ 30 -o RA14_TU_CL1-P0-gDNA.Tumor.sortedByCoord.bam - 
/usr/local/bin/gatk --java-options "-Xms30G" \
            MarkDuplicates \
            --INPUT RA14_TU_CL1-P0-gDNA.Normal.sortedByCoord.bam \
            --OUTPUT RA14_TU_CL1-P0-gDNA.Normal.duplicates_marked.bam \
            --METRICS_FILE RA14_TU_CL1-P0-gDNA.Normal.duplicate_metrics \
            --VALIDATION_STRINGENCY SILENT \
            --OPTICAL_DUPLICATE_PIXEL_DISTANCE 2500 \
            --ASSUME_SORT_ORDER "queryname" \
            --CREATE_MD5_FILE true 
/usr/local/bin/gatk --java-options "-Xms30G" \
            MarkDuplicates \
            --INPUT RA14_TU_CL1-P0-gDNA.Tumor.sortedByCoord.bam \
            --OUTPUT RA14_TU_CL1-P0-gDNA.Tumor.duplicates_marked.bam \
            --METRICS_FILE RA14_TU_CL1-P0-gDNA.Tumor.duplicate_metrics \
            --VALIDATION_STRINGENCY SILENT \
            --OPTICAL_DUPLICATE_PIXEL_DISTANCE 2500 \
            --ASSUME_SORT_ORDER "queryname" \
            --CREATE_MD5_FILE true 
/usr/local/bin/gatk SortSam\
            --INPUT RA14_TU_CL1-P0-gDNA.Normal.duplicates_marked.bam \
            --OUTPUT RA14_TU_CL1-P0-gDNA.Normal.duplicates_marked_sorted.bam \
            --SORT_ORDER "coordinate" \
            --CREATE_INDEX true \
            --CREATE_MD5_FILE false 
/usr/local/bin/gatk SetNmMdAndUqTags \
            --INPUT  RA14_TU_CL1-P0-gDNA.Normal.duplicates_marked_sorted.bam\
            --OUTPUT RA14_TU_CL1-P0-gDNA.Normal.duplicates_marked_sorted_fixed.bam \
            --CREATE_INDEX true \
            --CREATE_MD5_FILE true \
            --REFERENCE_SEQUENCE /cygene/database/GATK_resource_bundle/hg38/Homo_sapiens_assembly38.fasta 
/usr/local/bin/gatk SortSam\
            --INPUT RA14_TU_CL1-P0-gDNA.Tumor.duplicates_marked.bam \
            --OUTPUT RA14_TU_CL1-P0-gDNA.Tumor.duplicates_marked_sorted.bam \
            --SORT_ORDER "coordinate" \
            --CREATE_INDEX true \
            --CREATE_MD5_FILE false 
/usr/local/bin/gatk SetNmMdAndUqTags \
            --INPUT  RA14_TU_CL1-P0-gDNA.Tumor.duplicates_marked_sorted.bam \
            --OUTPUT RA14_TU_CL1-P0-gDNA.Tumor.duplicates_marked_sorted_fixed.bam \
            --CREATE_INDEX true \
            --CREATE_MD5_FILE true \
            --REFERENCE_SEQUENCE /cygene/database/GATK_resource_bundle/hg38/Homo_sapiens_assembly38.fasta 
/usr/local/bin/gatk --java-options "-Xms30G" \
            BaseRecalibrator \
            -R /cygene/database/GATK_resource_bundle/hg38/Homo_sapiens_assembly38.fasta \
            -I RA14_TU_CL1-P0-gDNA.Normal.duplicates_marked_sorted_fixed.bam \
            --use-original-qualities \
            -O RA14_TU_CL1-P0-gDNA.Normal.recal_data.csv \
            --known-sites /cygene/database/GATK_resource_bundle/hg38/dbsnp_138.hg38.vcf.gz \
            --known-sites /cygene/database/GATK_resource_bundle/hg38/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz \
            --known-sites /cygene/database/GATK_resource_bundle/hg38/hapmap_3.3.hg38.vcf.gz \
            --known-sites /cygene/database/GATK_resource_bundle/hg38/1000G_omni2.5.hg38.vcf.gz \
            --known-sites /cygene/database/GATK_resource_bundle/hg38/1000G_phase1.snps.high_confidence.hg38.vcf.gz
/usr/local/bin/gatk --java-options "-Xms30G" \
            BaseRecalibrator \
            -R /cygene/database/GATK_resource_bundle/hg38/Homo_sapiens_assembly38.fasta \
            -I RA14_TU_CL1-P0-gDNA.Tumor.duplicates_marked_sorted_fixed.bam \
            --use-original-qualities \
            -O RA14_TU_CL1-P0-gDNA.Tumor.recal_data.csv \
            --known-sites /cygene/database/GATK_resource_bundle/hg38/dbsnp_138.hg38.vcf.gz \
            --known-sites /cygene/database/GATK_resource_bundle/hg38/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz \
            --known-sites /cygene/database/GATK_resource_bundle/hg38/hapmap_3.3.hg38.vcf.gz \
            --known-sites /cygene/database/GATK_resource_bundle/hg38/1000G_omni2.5.hg38.vcf.gz \
            --known-sites /cygene/database/GATK_resource_bundle/hg38/1000G_phase1.snps.high_confidence.hg38.vcf.gz
/usr/local/bin/gatk --java-options "-Xms30G" \
            ApplyBQSR \
            -R /cygene/database/GATK_resource_bundle/hg38/Homo_sapiens_assembly38.fasta \
            -I  RA14_TU_CL1-P0-gDNA.Normal.duplicates_marked_sorted_fixed.bam \
            -O  RA14_TU_CL1-P0-gDNA.Normal.duplicates_marked_sorted_fixed.BQSR.bam \
            -bqsr RA14_TU_CL1-P0-gDNA.Normal.recal_data.csv \
            --create-output-bam-index \
            --static-quantized-quals 10 \
            --static-quantized-quals 20 \
            --static-quantized-quals 30 \
            --add-output-sam-program-record \
            --use-original-qualities \
            --create-output-bam-md5 true
/usr/local/bin/gatk --java-options "-Xms30G"  \
            ApplyBQSR \
            -R /cygene/database/GATK_resource_bundle/hg38/Homo_sapiens_assembly38.fasta \
            -I  RA14_TU_CL1-P0-gDNA.Tumor.duplicates_marked_sorted_fixed.bam  \
            -O  RA14_TU_CL1-P0-gDNA.Tumor.duplicates_marked_sorted_fixed.BQSR.bam \
            -bqsr RA14_TU_CL1-P0-gDNA.Tumor.recal_data.csv  \
            --static-quantized-quals 10 \
            --static-quantized-quals 20 \
            --static-quantized-quals 30 \
            --add-output-sam-program-record \
            --create-output-bam-index \
            --use-original-qualities \
            --create-output-bam-md5 true 
/usr/local/bin/gatk --java-options "-Xmx30G" \
            DepthOfCoverage \
            --input RA14_TU_CL1-P0-gDNA.Tumor.duplicates_marked_sorted_fixed.BQSR.bam \
            -L /cygene/work/00.test/pipeline/WES_cnv_somatic_pair_pipeline/database/whole_exome_illumina_hg38.targets.interval_list \
            -O RA14_TU_CL1-P0-gDNA.Tumor.bam.DepthOfCoverage.txt \
            -R /cygene/database/GATK_resource_bundle/hg38/Homo_sapiens_assembly38.fasta
/usr/local/bin/gatk --java-options "-Xmx30G" \
            DepthOfCoverage \
            --input RA14_TU_CL1-P0-gDNA.Normal.duplicates_marked_sorted_fixed.BQSR.bam  \
            -L /cygene/work/00.test/pipeline/WES_cnv_somatic_pair_pipeline/database/whole_exome_illumina_hg38.targets.interval_list \
            -O RA14_TU_CL1-P0-gDNA.Normal.bam.DepthOfCoverage.txt \
            -R /cygene/database/GATK_resource_bundle/hg38/Homo_sapiens_assembly38.fasta
Rscript /cygene/work/00.test/pipeline/WES_pipeline/bin/DepthOfCoverage.step1.deal.data.R RA14_TU_CL1-P0-gDNA.Tumor.bam.DepthOfCoverage.txt RA14_TU_CL1-P0-gDNA.Normal.bam.DepthOfCoverage.txt 
Rscript /cygene/work/00.test/pipeline/WES_pipeline/bin/DepthOfCoverage.step2.draw.plot.R RA14_TU_CL1-P0-gDNA.Tumor.coverage.depth.rate.xls RA14_TU_CL1-P0-gDNA.Normal.coverage.depth.rate.xls 
/usr/local/bin/gatk --java-options "-Xmx30G"  \
            Mutect2 \
            -R /cygene/database/GATK_resource_bundle/hg38/Homo_sapiens_assembly38.fasta \
            -I RA14_TU_CL1-P0-gDNA.Tumor.duplicates_marked_sorted_fixed.BQSR.bam \
            -I RA14_TU_CL1-P0-gDNA.Normal.duplicates_marked_sorted_fixed.BQSR.bam \
            -tumor Tumor \
            -normal Normal \
            -germline-resource /cygene/database/GATK_resource_bundle/Mutect2/af-only-gnomad.hg38.vcf.gz \
            -pon /cygene/database/GATK_resource_bundle/Mutect2/1000g_pon.hg38.vcf.gz   \
            -O RA14_TU_CL1-P0-gDNA.unfiltered.vcf \
            --f1r2-tar-gz RA14_TU_CL1-P0-gDNA.f1r2.tar.gz
/usr/local/bin/gatk --java-options "-Xmx30G"  \
            LearnReadOrientationModel \
            -I RA14_TU_CL1-P0-gDNA.f1r2.tar.gz \
            -O RA14_TU_CL1-P0-gDNA.read-orientation-model.tar.gz 
/usr/local/bin/gatk --java-options "-Xmx30G"  \
            GetPileupSummaries \
            -I RA14_TU_CL1-P0-gDNA.Tumor.duplicates_marked_sorted_fixed.BQSR.bam  \
            -V /cygene/database/GATK_resource_bundle/Mutect2/GetPileupSummaries/small_exac_common_3.hg38.vcf.gz \
            -L /cygene/database/GATK_resource_bundle/Mutect2/GetPileupSummaries/small_exac_common_3.hg38.vcf.gz \
            -O RA14_TU_CL1-P0-gDNA.Tumor.pileups.table 
/usr/local/bin/gatk --java-options "-Xmx30G"  \
            GetPileupSummaries \
            -I RA14_TU_CL1-P0-gDNA.Normal.duplicates_marked_sorted_fixed.BQSR.bam  \
            -V /cygene/database/GATK_resource_bundle/Mutect2/GetPileupSummaries/small_exac_common_3.hg38.vcf.gz \
            -L /cygene/database/GATK_resource_bundle/Mutect2/GetPileupSummaries/small_exac_common_3.hg38.vcf.gz \
            -O RA14_TU_CL1-P0-gDNA.Normal.pileups.table 
/usr/local/bin/gatk --java-options "-Xmx30G"  \
            CalculateContamination \
            -I RA14_TU_CL1-P0-gDNA.Tumor.pileups.table  \
            -matched RA14_TU_CL1-P0-gDNA.Normal.pileups.table \
            -O RA14_TU_CL1-P0-gDNA.Tumor.contamination.table \
            --tumor-segmentation RA14_TU_CL1-P0-gDNA.segments.table 
/usr/local/bin/gatk --java-options "-Xmx30G"  \
            FilterMutectCalls \
            -R /cygene/database/GATK_resource_bundle/hg38/Homo_sapiens_assembly38.fasta \
            -V RA14_TU_CL1-P0-gDNA.unfiltered.vcf \
            --contamination-table RA14_TU_CL1-P0-gDNA.Tumor.contamination.table \
            --tumor-segmentation RA14_TU_CL1-P0-gDNA.segments.table \
            --ob-priors RA14_TU_CL1-P0-gDNA.read-orientation-model.tar.gz \
            -stats RA14_TU_CL1-P0-gDNA.unfiltered.vcf.stats \
            --filtering-stats RA14_TU_CL1-P0-gDNA.filtering.stats\
            -O RA14_TU_CL1-P0-gDNA.filtered.vcf 
/usr/local/bin/gatk --java-options "-Xmx30G"  \
            Funcotator \
            --data-sources-path /cygene/database/GATK_resource_bundle/funcotator_dataSources.v1.7.20200521s \
            --ref-version hg38 \
            --output-file-format MAF \
            --reference /cygene/database/GATK_resource_bundle/hg38/Homo_sapiens_assembly38.fasta \
            --variant RA14_TU_CL1-P0-gDNA.filtered.vcf \
            --output RA14_TU_CL1-P0-gDNA.variants.funcotated.MAF.xls \
            --remove-filtered-variants true \
            --add-output-vcf-command-line false \
            --annotation-default normal_barcode:Normal \
            --annotation-default tumor_barcode:Tumor \
            --annotation-default Center:RootPath \
            --annotation-default Sequencer:Miseq 
grep -v "^#" RA14_TU_CL1-P0-gDNA.variants.funcotated.MAF.xls > RA14_TU_CL1-P0-gDNA.variants.funcotated.without.header.MAF.xls
python3 /cygene/work/00.test/pipeline/WES_pipeline/bin/extract_minigene.py /cygene/database/GATK_resource_bundle/hg38/gencode.v34.pc_transcripts.fa RA14_TU_CL1-P0-gDNA.variants.funcotated.without.header.MAF.xls RA14_TU_CL1-P0-gDNA.variants.funcotated.with.minigene.MAF.xls
less RA14_TU_CL1-P0-gDNA.variants.funcotated.with.minigene.MAF.xls | grep -v "Hugo_Symbol" |awk '{if ($5 == "MT") {print"chrM\t"$6-1"\t"$7}else{print$5"\t"$6-1"\t"$7}}' >RA14_TU_CL1-P0-gDNA.snp.checked.bed
/usr/local/bin/gatk --java-options "-Xmx30G" CollectReadCounts \
            -L /cygene/work/00.test/pipeline/GATK_somatic_CNV/database/hg38.filtered.interval_list \
            --input RA14_TU_CL1-P0-gDNA.Tumor.duplicates_marked_sorted_fixed.BQSR.bam \
            --reference /cygene/database/GATK_resource_bundle/hg38/Homo_sapiens_assembly38.fasta \
            --format HDF5 \
            --interval-merging-rule OVERLAPPING_ONLY \
            --output RA14_TU_CL1-P0-gDNA.Tumor.counts 
/usr/local/bin/gatk --java-options "-Xmx30G" CollectReadCounts \
            -L /cygene/work/00.test/pipeline/GATK_somatic_CNV/database/hg38.filtered.interval_list \
            --input RA14_TU_CL1-P0-gDNA.Normal.duplicates_marked_sorted_fixed.BQSR.bam \
            --reference /cygene/database/GATK_resource_bundle/hg38/Homo_sapiens_assembly38.fasta \
            --format HDF5 \
            --interval-merging-rule OVERLAPPING_ONLY \
            --output RA14_TU_CL1-P0-gDNA.Normal.counts 
/usr/local/bin/gatk --java-options "-Xmx30G" CreateReadCountPanelOfNormals \
            --input RA14_TU_CL1-P0-gDNA.Normal.counts \
            --minimum-interval-median-percentile 10.0 \
            --maximum-zeros-in-sample-percentage 5.0 \
            --maximum-zeros-in-interval-percentage 5.0 \
            --extreme-sample-median-percentile 2.5 \
            --do-impute-zeros true \
            --extreme-outlier-truncation-percentile 0.1 \
            --number-of-eigensamples 20 \
            --maximum-chunk-size 16777216 \
            --annotated-intervals /cygene/work/00.test/pipeline/GATK_somatic_CNV/database/hg38.annotated.tsv \
            --output RA14_TU_CL1-P0-gDNA.pon_entity_id.hdf5
/usr/local/bin/gatk --java-options "-Xmx30G" CollectAllelicCounts \
            -L /cygene/work/00.test/pipeline/GATK_somatic_CNV/database/hg38.filtered.interval_list \
            --input RA14_TU_CL1-P0-gDNA.Tumor.duplicates_marked_sorted_fixed.BQSR.bam \
            --reference /cygene/database/GATK_resource_bundle/hg38/Homo_sapiens_assembly38.fasta \
            --minimum-base-quality 20 \
            --output RA14_TU_CL1-P0-gDNA.Tumor.allelic_counts_file.txt 
/usr/local/bin/gatk --java-options "-Xmx30G" CollectAllelicCounts \
            -L /cygene/work/00.test/pipeline/GATK_somatic_CNV/database/hg38.filtered.interval_list \
            --input RA14_TU_CL1-P0-gDNA.Normal.duplicates_marked_sorted_fixed.BQSR.bam \
            --reference /cygene/database/GATK_resource_bundle/hg38/Homo_sapiens_assembly38.fasta \
            --minimum-base-quality 20 \
            --output RA14_TU_CL1-P0-gDNA.Normal.allelic_counts_file.txt 
/usr/local/bin/gatk --java-options "-Xms30G" DenoiseReadCounts \
            --INPUT RA14_TU_CL1-P0-gDNA.Tumor.read_counts_file.txt \
            --count-panel-of-normals RA14_TU_CL1-P0-gDNA.pon_entity_id.hdf5 \
            --number-of-eigensamples 20 \
            --standardized-copy-ratios RA14_TU_CL1-P0-gDNA.Tumor.standardizedCR.tsv \
            --denoised-copy-ratios RA14_TU_CL1-P0-gDNA.Tumor.denoisedCR.tsv 
/usr/local/bin/gatk --java-options "-Xms30G" DenoiseReadCounts \
            --INPUT RA14_TU_CL1-P0-gDNA.Normal.read_counts_file.txt \
            --count-panel-of-normals RA14_TU_CL1-P0-gDNA.pon_entity_id.hdf5 \
            --number-of-eigensamples 20 \
            --standardized-copy-ratios RA14_TU_CL1-P0-gDNA.Normal.standardizedCR.tsv \
            --denoised-copy-ratios RA14_TU_CL1-P0-gDNA.Normal.denoisedCR.tsv  
/usr/local/bin/gatk --java-options "-Xmx30G" ModelSegments \
            --denoised-copy-ratios RA14_TU_CL1-P0-gDNA.Tumor.denoisedCR.tsv \
            --allelic-counts RA14_TU_CL1-P0-gDNA.Tumor.allelic_counts_file.txt  \
            --normal-allelic-counts  RA14_TU_CL1-P0-gDNA.Normal.allelic_counts_file.txt \
            --minimum-total-allele-count-case 10 \
            --minimum-total-allele-count-normal 30 \
            --genotyping-homozygous-log-ratio-threshold "-10.0" \
            --genotyping-base-error-rate 0.05 \
            --maximum-number-of-segments-per-chromosome 1000 \
            --kernel-variance-copy-ratio  0.0 \
            --kernel-variance-allele-fraction 0.025 \
            --kernel-scaling-allele-fraction 1.0 \
            --kernel-approximation-dimension 100 \
            --window-size 256 \
            --number-of-changepoints-penalty-factor 1.0 \
            --minor-allele-fraction-prior-alpha 25.0 \
            --number-of-samples-copy-ratio 100 \
            --number-of-burn-in-samples-copy-ratio 50 \
            --number-of-samples-allele-fraction 100 \
            --number-of-burn-in-samples-allele-fraction 50 \
            --smoothing-credible-interval-threshold-copy-ratio 2.0\
            --smoothing-credible-interval-threshold-allele-fraction 2.0 \
            --maximum-number-of-smoothing-iterations 10 \
            --number-of-smoothing-iterations-per-fit 0 \
            --output ModelSegmentsTumor \
            --output-prefix RA14_TU_CL1-P0-gDNA.Tumor 
/usr/local/bin/gatk --java-options "-Xmx30G" ModelSegments \
            --denoised-copy-ratios RA14_TU_CL1-P0-gDNA.Normal.denoisedCR.tsv \
            --allelic-counts RA14_TU_CL1-P0-gDNA.Normal.allelic_counts_file.txt \
            --minimum-total-allele-count-case 5 \
            --minimum-total-allele-count-normal 30 \
            --genotyping-homozygous-log-ratio-threshold -10.0 \
            --genotyping-base-error-rate 0.05 \
            --maximum-number-of-segments-per-chromosome 1000 \
            --kernel-variance-copy-ratio 0.0 \
            --kernel-variance-allele-fraction 0.025 \
            --kernel-scaling-allele-fraction 1.0 \
            --kernel-approximation-dimension 100 \
            --window-size 256 \
            --number-of-changepoints-penalty-factor 1.0\
            --minor-allele-fraction-prior-alpha 25.0 \
            --number-of-samples-copy-ratio 100 \
            --number-of-burn-in-samples-copy-ratio 50 \
            --number-of-samples-allele-fraction 100 \
            --number-of-burn-in-samples-allele-fraction 50 \
            --smoothing-credible-interval-threshold-copy-ratio 2.0 \
            --smoothing-credible-interval-threshold-allele-fraction 2.0 \
            --maximum-number-of-smoothing-iterations 10 \
            --number-of-smoothing-iterations-per-fit 0 \
            --output ModelSegmentsNormal \
            --output-prefix RA14_TU_CL1-P0-gDNA.Normal 
