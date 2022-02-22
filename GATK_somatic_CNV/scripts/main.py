#!/usr/bin/python3.8
import os
import sys
import argparse
import configparser


def _argparse():
    parser = argparse.ArgumentParser(description="This is description")
    parser.add_argument('-c', '--config', action='store', dest='config',
                        default='config.txt', help="this is config file")
    # parser.add_argument('-l', '--list', action='store', dest='data', default='data.xls', help="this is data list file")
    return parser.parse_args()


def make_dir(*dir):
    for each in dir:
        if not os.path.exists(each):
            os.mkdir(each)


def main():
    """docstring for __main__"""
    parser = _argparse()
    cf = configparser.ConfigParser(allow_no_value=True)
    cf.read(parser.config)
    # check config file
    if not cf.has_section('Config'):
        os.exit("Error: your config file is not correct.")

    # read config file
    config_dict = {
        'gatk': cf.get('GATK', 'gatk'),
        'java_mem': cf.get('GATK', 'java_mem'),

        'known_dbSNP_vcf': cf.get('GATK', 'known_dbSNP_vcf'),
        'known_indel_sites_VCF': cf.get('GATK', 'known_indel_sites_VCF'),
        'known_hapmap_vcf': cf.get('GATK', 'known_hapmap_vcf'),
        'known_omni_vcf': cf.get('GATK', 'known_omni_vcf'),
        'known_1000G_phase1_snps_vcf': cf.get('GATK', 'known_1000G_phase1_snps_vcf'),

        # Mutect 2
        'af_only_gnomad_vcf': cf.get('GATK', 'af_only_gnomad_vcf'),
        'panel_of_normals_vcf': cf.get('GATK', 'panel_of_normals_vcf'),
        'variants_for_contamination': cf.get('GATK', 'variants_for_contamination'),

        # Funcotator
        'funcotator_dataSources': cf.get("GATK", "funcotator_dataSources"),
        'ref_version': cf.get("GATK", "ref_version"),

        # 'python' : cf.get('Config', 'python'),
        'bwa': cf.get('Config', 'bwa'),
        'bwa_threads': cf.get("Config", "bwa_threads"),
        'samtools': cf.get('Config', 'samtools'),
        'samtools_threads': cf.get('Config', 'samtools_threads'),
        'scripts_dir': cf.get('Config', 'scripts_dir'),

        # database
        'ref_fasta': cf.get('Config', 'ref_fasta'),
        'ref_fasta_dict': cf.get('Config', 'ref_fasta_dict'),
        'cds_fasta': cf.get('Config', 'cds_fasta'),

        'intervals': cf.get('Config', "intervals"),
        'annotated_intervals':cf.get('Config', "annotated_intervals"),
        'project_name': cf.get("Config", "project_name"),
        'sample_name': cf.get("Config", "sample_name"),
        'tumor_fastq': cf.get('Config', "tumor_fastq"),
        'normal_fastq': cf.get('Config', "normal_fastq"),
    }
    # make directories:
    project_dir = os.path.abspath(".") + '/' + config_dict['project_name']
    make_dir(project_dir)
    print("# Create work directory")

    # generate shell
    shell_name = project_dir + '/work.' + \
        config_dict['project_name'] + '.WES.sh'
    #shell_name = shell_dir + '/work.' + config_dict['project_name'] + '.sh'
    # only open a file so use try:finally to close.
    # rawdata_dict = deal_rawdata(parser.data, data_dir)

    with open(shell_name, "w") as f:
        # bwa.
        f.write("{bwa} mem -t {bwa_threads} -Y -H \"@HD\\tVN:1.5\\tGO:none\\tSO:coordinate\" -R \"@RG\\tID:{sample_name}_Normal\\tSM:Normal\\tLB:WES\\tPL:ILLumina\\tPU:HVW2MCCXX:6:none\" {ref_fasta} {normal_fastq} | {samtools} view -@ {samtools_threads} -buhS -t {ref_fasta}.fai - | {samtools} sort -@ {samtools_threads} -o {sample_name}.Normal.sortedByCoord.bam - \n".format(**config_dict))
        f.write("{bwa} mem -t {bwa_threads} -Y -H \"@HD\\tVN:1.5\\tGO:none\\tSO:coordinate\" -R \"@RG\\tID:{sample_name}_Tumor\\tSM:Tumor\\tLB:WES\\tPL:ILLumina\\tPU:HVW2MCCXX:6:none\" {ref_fasta} {tumor_fastq} | {samtools} view -@ {samtools_threads} -buhS -t {ref_fasta}.fai - | {samtools} sort -@ {samtools_threads} -o {sample_name}.Tumor.sortedByCoord.bam - \n".format(**config_dict))
        # MarkDuplicates
        f.write("""{gatk} --java-options \"-Xms{java_mem}G\" \\
            MarkDuplicates \\
            --INPUT {sample_name}.Normal.sortedByCoord.bam \\
            --OUTPUT {sample_name}.Normal.duplicates_marked.bam \\
            --METRICS_FILE {sample_name}.Normal.duplicate_metrics \\
            --VALIDATION_STRINGENCY SILENT \\
            --OPTICAL_DUPLICATE_PIXEL_DISTANCE 2500 \\
            --ASSUME_SORT_ORDER \"queryname\" \\
            --CREATE_MD5_FILE true \n""".format(**config_dict))
        f.write("""{gatk} --java-options \"-Xms{java_mem}G\" \\
            MarkDuplicates \\
            --INPUT {sample_name}.Tumor.sortedByCoord.bam \\
            --OUTPUT {sample_name}.Tumor.duplicates_marked.bam \\
            --METRICS_FILE {sample_name}.Tumor.duplicate_metrics \\
            --VALIDATION_STRINGENCY SILENT \\
            --OPTICAL_DUPLICATE_PIXEL_DISTANCE 2500 \\
            --ASSUME_SORT_ORDER \"queryname\" \\
            --CREATE_MD5_FILE true \n""".format(**config_dict))
        f.write("""{gatk} SortSam\\
            --INPUT {sample_name}.Normal.duplicates_marked.bam \\
            --OUTPUT {sample_name}.Normal.duplicates_marked_sorted.bam \\
            --SORT_ORDER "coordinate" \\
            --CREATE_INDEX true \\
            --CREATE_MD5_FILE false \n""".format(**config_dict))
        f.write("""{gatk} SetNmMdAndUqTags \\
            --INPUT  {sample_name}.Normal.duplicates_marked_sorted.bam\\
            --OUTPUT {sample_name}.Normal.duplicates_marked_sorted_fixed.bam \\
            --CREATE_INDEX true \\
            --CREATE_MD5_FILE true \\
            --REFERENCE_SEQUENCE {ref_fasta} \n""".format(**config_dict))
        f.write("""{gatk} SortSam\\
            --INPUT {sample_name}.Tumor.duplicates_marked.bam \\
            --OUTPUT {sample_name}.Tumor.duplicates_marked_sorted.bam \\
            --SORT_ORDER "coordinate" \\
            --CREATE_INDEX true \\
            --CREATE_MD5_FILE false \n""".format(**config_dict))
        f.write("""{gatk} SetNmMdAndUqTags \\
            --INPUT  {sample_name}.Tumor.duplicates_marked_sorted.bam \\
            --OUTPUT {sample_name}.Tumor.duplicates_marked_sorted_fixed.bam \\
            --CREATE_INDEX true \\
            --CREATE_MD5_FILE true \\
            --REFERENCE_SEQUENCE {ref_fasta} \n""".format(**config_dict))
        # # BaseRecalibrator
        f.write("""{gatk} --java-options \"-Xms{java_mem}G\" \\
            BaseRecalibrator \\
            -R {ref_fasta} \\
            -I {sample_name}.Normal.duplicates_marked_sorted_fixed.bam \\
            --use-original-qualities \\
            -O {sample_name}.Normal.recal_data.csv \\
            --known-sites {known_dbSNP_vcf} \\
            --known-sites {known_indel_sites_VCF} \\
            --known-sites {known_hapmap_vcf} \\
            --known-sites {known_omni_vcf} \\
            --known-sites {known_1000G_phase1_snps_vcf}\n""".format(**config_dict))
        f.write("""{gatk} --java-options \"-Xms{java_mem}G\" \\
            BaseRecalibrator \\
            -R {ref_fasta} \\
            -I {sample_name}.Tumor.duplicates_marked_sorted_fixed.bam \\
            --use-original-qualities \\
            -O {sample_name}.Tumor.recal_data.csv \\
            --known-sites {known_dbSNP_vcf} \\
            --known-sites {known_indel_sites_VCF} \\
            --known-sites {known_hapmap_vcf} \\
            --known-sites {known_omni_vcf} \\
            --known-sites {known_1000G_phase1_snps_vcf}\n""".format(**config_dict))
        # # ApplylyBQSR
        f.write("""{gatk} --java-options \"-Xms{java_mem}G\" \\
            ApplyBQSR \\
            -R {ref_fasta} \\
            -I  {sample_name}.Normal.duplicates_marked_sorted_fixed.bam \\
            -O  {sample_name}.Normal.duplicates_marked_sorted_fixed.BQSR.bam \\
            -bqsr {sample_name}.Normal.recal_data.csv \\
            --create-output-bam-index \\
            --static-quantized-quals 10 \\
            --static-quantized-quals 20 \\
            --static-quantized-quals 30 \\
            --add-output-sam-program-record \\
            --use-original-qualities \\
            --create-output-bam-md5 true\n""".format(**config_dict))
        f.write("""{gatk} --java-options \"-Xms{java_mem}G\"  \\
            ApplyBQSR \\
            -R {ref_fasta} \\
            -I  {sample_name}.Tumor.duplicates_marked_sorted_fixed.bam  \\
            -O  {sample_name}.Tumor.duplicates_marked_sorted_fixed.BQSR.bam \\
            -bqsr {sample_name}.Tumor.recal_data.csv  \\
            --static-quantized-quals 10 \\
            --static-quantized-quals 20 \\
            --static-quantized-quals 30 \\
            --add-output-sam-program-record \\
            --create-output-bam-index \\
            --use-original-qualities \\
            --create-output-bam-md5 true \n""".format(**config_dict))
        # statistic bam
        f.write("""{gatk} --java-options \"-Xmx{java_mem}G\" \\
            DepthOfCoverage \\
            --input {sample_name}.Tumor.duplicates_marked_sorted_fixed.BQSR.bam \\
            -L /cygene/work/00.test/pipeline/WES_cnv_somatic_pair_pipeline/database/whole_exome_illumina_hg38.targets.interval_list \\
            -O {sample_name}.Tumor.bam.DepthOfCoverage.txt \\
            -R {ref_fasta}\n""".format(**config_dict))
        f.write("""{gatk} --java-options \"-Xmx{java_mem}G\" \\
            DepthOfCoverage \\
            --input {sample_name}.Normal.duplicates_marked_sorted_fixed.BQSR.bam  \\
            -L /cygene/work/00.test/pipeline/WES_cnv_somatic_pair_pipeline/database/whole_exome_illumina_hg38.targets.interval_list \\
            -O {sample_name}.Normal.bam.DepthOfCoverage.txt \\
            -R {ref_fasta}\n""".format(**config_dict))
        f.write("""Rscript {scripts_dir}/DepthOfCoverage.step1.deal.data.R {sample_name}.Tumor.bam.DepthOfCoverage.txt {sample_name}.Normal.bam.DepthOfCoverage.txt \n""".format(**config_dict))
        f.write("""Rscript {scripts_dir}/DepthOfCoverage.step2.draw.plot.R {sample_name}.Tumor.coverage.depth.rate.xls {sample_name}.Normal.coverage.depth.rate.xls \n""".format(**config_dict))

        # # Mutect2
        f.write("""{gatk} --java-options \"-Xmx{java_mem}G\"  \\
            Mutect2 \\
            -R {ref_fasta} \\
            -I {sample_name}.Tumor.duplicates_marked_sorted_fixed.BQSR.bam \\
            -I {sample_name}.Normal.duplicates_marked_sorted_fixed.BQSR.bam \\
            -tumor Tumor \\
            -normal Normal \\
            -germline-resource {af_only_gnomad_vcf} \\
            -pon {panel_of_normals_vcf}   \\
            -O {sample_name}.unfiltered.vcf \\
            --f1r2-tar-gz {sample_name}.f1r2.tar.gz\n""".format(**config_dict))
        f.write("""{gatk} --java-options \"-Xmx{java_mem}G\"  \\
            LearnReadOrientationModel \\
            -I {sample_name}.f1r2.tar.gz \\
            -O {sample_name}.read-orientation-model.tar.gz \n""".format(**config_dict))
        # contamination.
        f.write("""{gatk} --java-options \"-Xmx{java_mem}G\"  \\
            GetPileupSummaries \\
            -I {sample_name}.Tumor.duplicates_marked_sorted_fixed.BQSR.bam  \\
            -V {variants_for_contamination} \\
            -L {variants_for_contamination} \\
            -O {sample_name}.Tumor.pileups.table \n""".format(**config_dict))
        f.write("""{gatk} --java-options \"-Xmx{java_mem}G\"  \\
            GetPileupSummaries \\
            -I {sample_name}.Normal.duplicates_marked_sorted_fixed.BQSR.bam  \\
            -V {variants_for_contamination} \\
            -L {variants_for_contamination} \\
            -O {sample_name}.Normal.pileups.table \n""".format(**config_dict))
        f.write("""{gatk} --java-options \"-Xmx{java_mem}G\"  \\
            CalculateContamination \\
            -I {sample_name}.Tumor.pileups.table  \\
            -matched {sample_name}.Normal.pileups.table \\
            -O {sample_name}.Tumor.contamination.table \\
            --tumor-segmentation {sample_name}.segments.table \n""".format(**config_dict))

        # filter mutation
        f.write("""{gatk} --java-options \"-Xmx{java_mem}G\"  \\
            FilterMutectCalls \\
            -R {ref_fasta} \\
            -V {sample_name}.unfiltered.vcf \\
            --contamination-table {sample_name}.Tumor.contamination.table \\
            --tumor-segmentation {sample_name}.segments.table \\
            --ob-priors {sample_name}.read-orientation-model.tar.gz \\
            -stats {sample_name}.unfiltered.vcf.stats \\
            --filtering-stats {sample_name}.filtering.stats\\
            -O {sample_name}.filtered.vcf \n""".format(**config_dict))

        # Function anatation.
        f.write("""{gatk} --java-options \"-Xmx{java_mem}G\"  \\
            Funcotator \\
            --data-sources-path {funcotator_dataSources} \\
            --ref-version {ref_version} \\
            --output-file-format MAF \\
            --reference {ref_fasta} \\
            --variant {sample_name}.filtered.vcf \\
            --output {sample_name}.variants.funcotated.MAF.xls \\
            --remove-filtered-variants true \\
            --add-output-vcf-command-line false \\
            --annotation-default normal_barcode:Normal \\
            --annotation-default tumor_barcode:Tumor \\
            --annotation-default Center:RootPath \\
            --annotation-default Sequencer:Miseq \n""".format(**config_dict))
        f.write(
            """grep -v \"^#\" {sample_name}.variants.funcotated.MAF.xls > {sample_name}.variants.funcotated.without.header.MAF.xls\n""".format(**config_dict))
        f.write("""python3 {scripts_dir}/extract_minigene.py {cds_fasta} {sample_name}.variants.funcotated.without.header.MAF.xls {sample_name}.variants.funcotated.with.minigene.MAF.xls\n""".format(**config_dict))
        #f.write("""less {sample_name}.variants.funcotated.with.minigene.MAF.xls | grep -v "Hugo_Symbol" |awk '{{print$5"\\t"$6-1"\\t"$7}}' > {sample_name}.snp.checked.bed\n""".format(**config_dict))
        f.write("""less {sample_name}.variants.funcotated.with.minigene.MAF.xls | grep -v "Hugo_Symbol" |awk '{{if ($5 == "MT") {{print"chrM\\t"$6-1"\\t"$7}}else{{print$5"\\t"$6-1"\\t"$7}}}}' >{sample_name}.snp.checked.bed\n""".format(**config_dict))

        ###########################################################
        # For CNVs
        ###########################################################

        # CNVTasks.PreprocessIntervals
        # AnnotateIntervals
        # FilterIntervals
        # CollectCountsTumor
        f.write("""{gatk} --java-options "-Xmx30G" CollectReadCounts \\
            -L {intervals} \\
            --input {sample_name}.Tumor.duplicates_marked_sorted_fixed.BQSR.bam \\
            --reference {ref_fasta} \\
            --format HDF5 \\
            --interval-merging-rule OVERLAPPING_ONLY \\
            --output {sample_name}.Tumor.counts \n""".format(**config_dict))
        # CollectCountsNormal
        f.write("""{gatk} --java-options "-Xmx30G" CollectReadCounts \\
            -L {intervals} \\
            --input {sample_name}.Normal.duplicates_marked_sorted_fixed.BQSR.bam \\
            --reference {ref_fasta} \\
            --format HDF5 \\
            --interval-merging-rule OVERLAPPING_ONLY \\
            --output {sample_name}.Normal.counts \n""".format(**config_dict))

        # CreateReadCountPanelOfNormals
        f.write("""{gatk} --java-options "-Xmx{java_mem}G" CreateReadCountPanelOfNormals \\
            --input {sample_name}.Normal.counts \\
            --minimum-interval-median-percentile 10.0 \\
            --maximum-zeros-in-sample-percentage 5.0 \\
            --maximum-zeros-in-interval-percentage 5.0 \\
            --extreme-sample-median-percentile 2.5 \\
            --do-impute-zeros true \\
            --extreme-outlier-truncation-percentile 0.1 \\
            --number-of-eigensamples 20 \\
            --maximum-chunk-size 16777216 \\
            --annotated-intervals {annotated_intervals} \\
            --output {sample_name}.pon_entity_id.hdf5\n""".format(**config_dict))

        # CollectAllelicCountsTumor
        f.write("""{gatk} --java-options "-Xmx{java_mem}G" CollectAllelicCounts \\
            -L {intervals} \\
            --input {sample_name}.Tumor.duplicates_marked_sorted_fixed.BQSR.bam \\
            --reference {ref_fasta} \\
            --minimum-base-quality 20 \\
            --output {sample_name}.Tumor.allelic_counts_file.txt \n""".format(**config_dict))
        # CollectAllelicCountsNormal
        f.write("""{gatk} --java-options "-Xmx{java_mem}G" CollectAllelicCounts \\
            -L {intervals} \\
            --input {sample_name}.Normal.duplicates_marked_sorted_fixed.BQSR.bam \\
            --reference {ref_fasta} \\
            --minimum-base-quality 20 \\
            --output {sample_name}.Normal.allelic_counts_file.txt \n""".format(**config_dict))

        # DenoiseReadCountsTumor
        f.write("""{gatk} --java-options "-Xms{java_mem}G" DenoiseReadCounts \\
            --INPUT {sample_name}.Tumor.read_counts_file.txt \\
            --count-panel-of-normals {sample_name}.pon_entity_id.hdf5 \\
            --number-of-eigensamples 20 \\
            --standardized-copy-ratios {sample_name}.Tumor.standardizedCR.tsv \\
            --denoised-copy-ratios {sample_name}.Tumor.denoisedCR.tsv \n""".format(**config_dict))
        # DenoiseReadCountsNormal
        f.write("""{gatk} --java-options "-Xms{java_mem}G" DenoiseReadCounts \\
            --INPUT {Normal.read_counts} \\
            --count-panel-of-normals {sample_name}.pon_entity_id.hdf5 \\
            --number-of-eigensamples 20 \\
            --standardized-copy-ratios {sample_name}.Normal.standardizedCR.tsv \\
            --denoised-copy-ratios {sample_name}.Normal.denoisedCR.tsv  \n""".format(**config_dict))
        
        # ModelSegmentsTumor
        f.write("""{gatk} --java-options "-Xmx30G" ModelSegments \\
            --denoised-copy-ratios {sample_name}.Tumor.denoisedCR.tsv \\
            --allelic-counts {sample_name}.Tumor.allelic_counts_file.txt  \\
            --normal-allelic-counts  {sample_name}.Normal.allelic_counts_file.txt \\
            --minimum-total-allele-count-case 10 \\
            --minimum-total-allele-count-normal 30 \\
            --genotyping-homozygous-log-ratio-threshold "-10.0" \\
            --genotyping-base-error-rate 0.05 \\
            --maximum-number-of-segments-per-chromosome 1000 \\
            --kernel-variance-copy-ratio  0.0 \\
            --kernel-variance-allele-fraction 0.025 \\
            --kernel-scaling-allele-fraction 1.0 \\
            --kernel-approximation-dimension 100 \\
            --window-size 256 \\
            --number-of-changepoints-penalty-factor 1.0 \\
            --minor-allele-fraction-prior-alpha 25.0 \\
            --number-of-samples-copy-ratio 100 \\
            --number-of-burn-in-samples-copy-ratio 50 \\
            --number-of-samples-allele-fraction 100 \\
            --number-of-burn-in-samples-allele-fraction 50 \\
            --smoothing-credible-interval-threshold-copy-ratio 2.0\\
            --smoothing-credible-interval-threshold-allele-fraction 2.0 \\
            --maximum-number-of-smoothing-iterations 10 \\
            --number-of-smoothing-iterations-per-fit 0 \\
            --output ModelSegmentsTumor \\
            --output-prefix {sample_name}.Tumor \n""".format(**config_dict))
        # ModelSegmentsNormal
        f.write("""{gatk} --java-options "-Xmx30G" \\
            ModelSegments \\
            --denoised-copy-ratios {sample_name}.Normal.denoisedCR.tsv \\
            --allelic-counts {sample_name}.Normal.allelic_counts_file.txt \\
            --minimum-total-allele-count-case 5 \\
            --minimum-total-allele-count-normal 30 \\
            --genotyping-homozygous-log-ratio-threshold -10.0 \\
            --genotyping-base-error-rate 0.05 \\
            --maximum-number-of-segments-per-chromosome 1000 \\
            --kernel-variance-copy-ratio 0.0 \\
            --kernel-variance-allele-fraction 0.025 \\
            --kernel-scaling-allele-fraction 1.0 \\
            --kernel-approximation-dimension 100 \\
            --window-size 256 \\
            --number-of-changepoints-penalty-factor 1.0\\
            --minor-allele-fraction-prior-alpha 25.0 \\
            --number-of-samples-copy-ratio 100 \\
            --number-of-burn-in-samples-copy-ratio 50 \\
            --number-of-samples-allele-fraction 100 \\
            --number-of-burn-in-samples-allele-fraction 50 \\
            --smoothing-credible-interval-threshold-copy-ratio 2.0 \\
            --smoothing-credible-interval-threshold-allele-fraction 2.0 \\
            --maximum-number-of-smoothing-iterations 10 \\
            --number-of-smoothing-iterations-per-fit 0 \\
            --output ModelSegmentsNormal \\
            --output-prefix {sample_name}.Normal \n""".format(**config_dict))
        
        # CallCopyRatioSegmentsTumor
        f.write("""{gatk} --java-options "-Xmx30G" CallCopyRatioSegments \\
            --input {copy_ratio_segments} \\
            --neutral-segment-copy-ratio-lower-bound 0.9\\
            --neutral-segment-copy-ratio-upper-bound 1.1 \\
            --outlier-neutral-segment-copy-ratio-z-score-threshold 2.0 \\
            --calling-copy-ratio-z-score-threshold 2.0 \\
            --output {sample_name}.called.seg \n""".format(**config_dict))
        # CallCopyRatioSegmentsNormal
        f.write("""{gatk} --java-options "-Xmx30G" CallCopyRatioSegments \\
            --input {copy_ratio_segments} \\
            --neutral-segment-copy-ratio-lower-bound 0.9 \
            --neutral-segment-copy-ratio-upper-bound 1.1 \
            --outlier-neutral-segment-copy-ratio-z-score-threshold 2.0 \
            --calling-copy-ratio-z-score-threshold 2.0 \
            --output {sample_name}.called.seg \n""".format(**config_dict))
        # output {
        #        File called_copy_ratio_segments = "{sample_name}.called.seg"
        #        File called_copy_ratio_legacy_segments = "{sample_name}.called.igv.seg"
        #     }
        # PlotDenoisedCopyRatiosTumor
        f.write("""{gatk} --java-options "-Xmx30G" PlotDenoisedCopyRatios \
            --standardized-copy-ratios {standardized_copy_ratios} \
            --denoised-copy-ratios {denoised_copy_ratios} \
            --sequence-dictionary {ref_fasta_dict} \
            --minimum-contig-length 1000000 \
            --output {output_dir_} \
            --output-prefix {sample_name} \n""".format(**config_dict))
        # PlotDenoisedCopyRatiosNormal
        f.write("""{gatk} --java-options "-Xmx{java_mem}G" PlotDenoisedCopyRatios \
            --standardized-copy-ratios {standardized_copy_ratios} \
            --denoised-copy-ratios {denoised_copy_ratios} \
            --sequence-dictionary {ref_fasta_dict} \
            --minimum-contig-length 1000000 \
            --output {output_dir_} \
            --output-prefix {sample_name} \n""".format(**config_dict))
        # PlotModeledSegmentsTumor
        f.write("""{gatk} --java-options "-Xmx{java_mem}G" PlotModeledSegments \
            --denoised-copy-ratios {denoised_copy_ratios} \
            --allelic-counts {het_allelic_counts} \
            --segments {modeled_segments} \
            --sequence-dictionary {ref_fasta_dict} \
            --minimum-contig-length 1000000 \
            --output {output_dir_} \
            --output-prefix {sample_name} \n""".format(**config_dict))
        # PlotModeledSegmentsNormal
        f.write("""{gatk} --java-options "-Xmx{java_mem}G" PlotModeledSegments \
            --denoised-copy-ratios {denoised_copy_ratios} \
            --allelic-counts {het_allelic_counts} \
            --segments {modeled_segments} \
            --sequence-dictionary {ref_fasta_dict} \
            --minimum-contig-length 1000000 \
            --output {output_dir_} \
            --output-prefix {sample_name} \n""".format(**config_dict))
        # CNVFuncotateSegments
        f.write("""{gatk} --java-options "-Xmx{java_mem}G" FuncotateSegments \
             --data-sources-path $DATA_SOURCES_FOLDER \
             --ref-version {funcotator_ref_version} \
             --output-file-format SEG \
             -R {ref_fasta} \
             --segments {input_seg_file} \
             -O {basename_input_seg_file}.funcotated.tsv \
             {interval_list_arg} " interval_list} \
             {"--transcript-selection-mode " + transcript_selection_mode} \
             {transcript_selection_arg}" sep=" --transcript-list " transcript_selection_list} \
             {annotation_def_arg}" sep=" --annotation-default " annotation_defaults} \
             {annotation_over_arg}" sep=" --annotation-override " annotation_overrides} \
             {excluded_fields_args}" sep=" --exclude-field " funcotator_excluded_fields} \
             {extra_args_arg}""".format(**config_dict))


print("all finished!")

if __name__ == '__main__':
    main()
