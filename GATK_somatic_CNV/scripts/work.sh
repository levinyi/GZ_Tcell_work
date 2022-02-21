# PON

# CNVTasks.PreprocessIntervals 
gatk --java-options "-Xmx30G" PreprocessIntervals \
            --sequence-dictionary ${ref_fasta_dict} \
            --reference ${ref_fasta} \
            --padding  250 \
            --bin-length 1000 \
            --interval-merging-rule OVERLAPPING_ONLY \
            --output ${base_filename}.preprocessed.interval_list
########### AnnotateIntervals
gatk --java-options "-Xmx${command_mem_mb}m" AnnotateIntervals \
            -L ${intervals} \
            --reference ${ref_fasta} \
            # ${"--mappability-track " + mappability_track_bed} \
            # ${"--segmental-duplication-track " + segmental_duplication_track_bed} \
            --feature-query-lookahead 1000000 \
            --interval-merging-rule OVERLAPPING_ONLY \
            --output ${base_filename}.annotated.tsv
########## FilterIntervals
gatk --java-options "-Xmx30G" FilterIntervals \
            -L ${intervals} \
            # ${"-XL " + blacklist_intervals} \
            --annotated-intervals ${base_filename}.annotated.tsv \
            # ${if defined(read_count_files) then "--input " else ""} ${sep=" --input " read_count_files} \
            --minimum-gc-content 0.1 \
            --maximum-gc-content 0.9 \
            --minimum-mappability 0.9 \
            --maximum-mappability 1.0 \
            --minimum-segmental-duplication-content 0.0 \
            --maximum-segmental-duplication-content 0.5 \
            --low-count-filter-count-threshold 5 \
            --low-count-filter-percentage-of-samples 90 \
            --extreme-count-filter-minimum-percentile 1.0 \
            --extreme-count-filter-maximum-percentile 99.0 \
            --extreme-count-filter-percentage-of-samples 90.0 \
            --interval-merging-rule OVERLAPPING_ONLY \
            --output ${base_filename}.filtered.interval_list
# CollectCountsTumor
gatk --java-options "-Xmx30G" CollectReadCounts \
            -L ${base_filename}.filtered.interval_list \
            --input ${Tumor.bam} \
            --reference ${ref_fasta} \
            --format HDF5 \
            --interval-merging-rule OVERLAPPING_ONLY \
            --output ${counts_filename}
# CollectCountsNormal
gatk --java-options "-Xmx30G" CollectReadCounts \
            -L ${base_filename}.filtered.interval_list \
            --input ${Normal.bam} \
            --reference ${ref_fasta} \
            --format HDF5 \
            --interval-merging-rule OVERLAPPING_ONLY \
            --output ${counts_filename}
# CollectAllelicCountsTumor
gatk --java-options "-Xmx30G" CollectAllelicCounts \
            -L ${common_sites} \ ??
            --input ${Tumor.bam} \
            --reference ${ref_fasta} \
            --minimum-base-quality 20 \
            --output ${Tumor.allelic_counts_filename}
# CollectAllelicCountsNormal
gatk --java-options "-Xmx30G" CollectAllelicCounts \
            -L ${common_sites} \  ??
            --input ${Normal.BQSR.bam} \
            --reference ${ref_fasta} \
            --minimum-base-quality 20 \
            --output ${allelic_counts_filename}

# DenoiseReadCountsTumor
gatk --java-options "-Xms30G" DenoiseReadCounts \
            --INPUT ${Tumor.read_counts} \ # from CollectReadCounts result
	    --count-panel-of-normals ${read_count_pon} \ 1000G pon file?
            # ${"--number-of-eigensamples " + number_of_eigensamples} \
            --standardized-copy-ratios ${entity_id}.standardizedCR.tsv \ # output
            --denoised-copy-ratios ${entity_id}.denoisedCR.tsv # output
# DenoiseReadCountsNormal
gatk --java-options "-Xms30G" DenoiseReadCounts \
            --INPUT ${Normal.read_counts} \
	    --count-panel-of-normals ${read_count_pon} \
            # ${"--number-of-eigensamples " + number_of_eigensamples} \
            --standardized-copy-ratios ${entity_id}.standardizedCR.tsv \ # output
            --denoised-copy-ratios ${entity_id}.denoisedCR.tsv # output
# ModelSegmentsTumor
gatk --java-options "-Xmx30G" ModelSegments \
	--denoised-copy-ratios ${denoised_copy_ratios} \
        --allelic-counts ${Tumor.allelic_counts_filename}  \
        --normal-allelic-counts  ${Normal.allelic_counts_filename} \
        --minimum-total-allele-count-case ${min_total_allele_count_} \??
        --minimum-total-allele-count-normal 30 \
        --genotyping-homozygous-log-ratio-threshold "-10.0" \
        --genotyping-base-error-rate 0.05 \
        --maximum-number-of-segments-per-chromosome 1000 \
        --kernel-variance-copy-ratio  0.0 \
        --kernel-variance-allele-fraction 0.025 \
        --kernel-scaling-allele-fraction 1.0 \
        --kernel-approximation-dimension 100 \
        --window-size ${window_sizes} \??
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
        --output ${output_dir_} \
        --output-prefix ${entity_id}
# ModelSegmentsNormal
gatk --java-options "-Xmx30G" \
	ModelSegments \
	--denoised-copy-ratios ${denoised_copy_ratios} \
        --allelic-counts ${allelic_counts} \
        ${normal_allelic_counts} \
        --minimum-total-allele-count-case ${min_total_allele_count_} \
        --minimum-total-allele-count-normal 30 \
        --genotyping-homozygous-log-ratio-threshold -10.0 \
        --genotyping-base-error-rate 0.05 \
        --maximum-number-of-segments-per-chromosome 1000 \
        --kernel-variance-copy-ratio 0.0 \
        --kernel-variance-allele-fraction 0.025 \
        --kernel-scaling-allele-fraction 1.0 \
        --kernel-approximation-dimension 100 \
        --window-size ${window_sizes} \
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
        --output ${output_dir_} \
        --output-prefix ${entity_id}
# CallCopyRatioSegmentsTumor
gatk --java-options "-Xmx30G" CallCopyRatioSegments \
            --input ${copy_ratio_segments} \
            --neutral-segment-copy-ratio-lower-bound 0.9\
            --neutral-segment-copy-ratio-upper-bound 1.1 \
            --outlier-neutral-segment-copy-ratio-z-score-threshold 2.0 \
            --calling-copy-ratio-z-score-threshold 2.0 \
            --output ${entity_id}.called.seg
# CallCopyRatioSegmentsNormal
gatk --java-options "-Xmx30G" CallCopyRatioSegments \
            --input ${copy_ratio_segments} \
            --neutral-segment-copy-ratio-lower-bound 0.9 \
            --neutral-segment-copy-ratio-upper-bound 1.1 \
            --outlier-neutral-segment-copy-ratio-z-score-threshold 2.0 \
            --calling-copy-ratio-z-score-threshold 2.0 \
            --output ${entity_id}.called.seg
# output {
#        File called_copy_ratio_segments = "${entity_id}.called.seg"
#        File called_copy_ratio_legacy_segments = "${entity_id}.called.igv.seg"
#     }
# PlotDenoisedCopyRatiosTumor
gatk --java-options "-Xmx30G" PlotDenoisedCopyRatios \
            --standardized-copy-ratios ${standardized_copy_ratios} \
            --denoised-copy-ratios ${denoised_copy_ratios} \
            --sequence-dictionary ${ref_fasta_dict} \
            --minimum-contig-length 1000000 \
            --output ${output_dir_} \
            --output-prefix ${entity_id}
# PlotDenoisedCopyRatiosNormal
gatk --java-options "-Xmx${command_mem_mb}m" PlotDenoisedCopyRatios \
            --standardized-copy-ratios ${standardized_copy_ratios} \
            --denoised-copy-ratios ${denoised_copy_ratios} \
            --sequence-dictionary ${ref_fasta_dict} \
            --minimum-contig-length 1000000 \
            --output ${output_dir_} \
            --output-prefix ${entity_id}
# PlotModeledSegmentsTumor
gatk --java-options "-Xmx${command_mem_mb}m" PlotModeledSegments \
            --denoised-copy-ratios ${denoised_copy_ratios} \
            --allelic-counts ${het_allelic_counts} \
            --segments ${modeled_segments} \
            --sequence-dictionary ${ref_fasta_dict} \
            --minimum-contig-length 1000000 \
            --output ${output_dir_} \
            --output-prefix ${entity_id}
# PlotModeledSegmentsNormal
gatk --java-options "-Xmx${command_mem_mb}m" PlotModeledSegments \
            --denoised-copy-ratios ${denoised_copy_ratios} \
            --allelic-counts ${het_allelic_counts} \
            --segments ${modeled_segments} \
            --sequence-dictionary ${ref_fasta_dict} \
            --minimum-contig-length 1000000 \
            --output ${output_dir_} \
            --output-prefix ${entity_id}
# CNVFuncotateSegments
gatk --java-options "-Xmx${command_mem_mb}m" FuncotateSegments \
             --data-sources-path $DATA_SOURCES_FOLDER \
             --ref-version ${funcotator_ref_version} \
             --output-file-format SEG \
             -R ${ref_fasta} \
             --segments ${input_seg_file} \
             -O ${basename_input_seg_file}.funcotated.tsv \
             ${interval_list_arg} " interval_list} \
             ${"--transcript-selection-mode " + transcript_selection_mode} \
             ${transcript_selection_arg}" sep=" --transcript-list " transcript_selection_list} \
             ${annotation_def_arg}" sep=" --annotation-default " annotation_defaults} \
             ${annotation_over_arg}" sep=" --annotation-override " annotation_overrides} \
             ${excluded_fields_args}" sep=" --exclude-field " funcotator_excluded_fields} \
             ${extra_args_arg}
