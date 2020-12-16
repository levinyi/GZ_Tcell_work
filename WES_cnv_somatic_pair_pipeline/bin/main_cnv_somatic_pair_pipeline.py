#!/usr/bin/python3.8
import os
import sys
import argparse
import configparser


def _argparse():
    parser = argparse.ArgumentParser(description="This is description")
    parser.add_argument('-c', '--config', action='store', dest='config', default='config.txt', help="this is config file")
    # parser.add_argument('-l', '--list', action='store', dest='data', default='data.xls', help="this is data list file")
    return parser.parse_args()


def make_dir(*dir):
    for each in dir:
        if not os.path.exists(each):
            os.mkdir( each )



def main():
    """docstring for __main__"""
    parser = _argparse()
    cf = configparser.ConfigParser(allow_no_value=True)
    cf.read(parser.config)
    # excel_file = parser.data
    # excel_file = os.path.abspath(excel_file)
    # check config file
    if not cf.has_section('Config'):
        os.exit("Error: your config file is not correct.")
    
    # read config file
    config_dict = {
        'gatk' : cf.get('GATK', 'gatk'),
        'java_mem': cf.get('GATK','java_mem'),

        # Funcotator
        'funcotator_dataSources' : cf.get("GATK", "funcotator_dataSources"),
        'ref_version' : cf.get("GATK", "ref_version"),
        
        # 'python' : cf.get('Config', 'python'),
        'scripts_dir': cf.get('Config', 'scripts_dir'),
        
        # database or files
        'ref_fasta': cf.get('Config', 'ref_fasta'),
        'common_sites': cf.get('GATK', 'common_sites'),
        'intervals'   : cf.get('GATK', "intervals"),
        'read_count_pon': cf.get('GATK', "read_count_pon"),
        # 'blacklist_intervals': cf.get('Config','blacklist_intervals'),
        'project_name' : cf.get("Config", "project_name"),
        'sample_name'  : cf.get("Config", "sample_name"),
        'tumor_bam'  : cf.get('Config', "tumor_bam"),
        'normal_bam' : cf.get('Config', "normal_bam"),
    }
    #make directories:
    project_dir = os.path.abspath(".") + '/' + config_dict['project_name']
    make_dir(project_dir)
    print("# Create work directory")

    # generate shell
    shell_name = project_dir + '/work.' + config_dict['project_name'] + '.WES.sh'
    #shell_name = shell_dir + '/work.' + config_dict['project_name'] + '.sh'
    # only open a file so use try:finally to close.
    # rawdata_dict = deal_rawdata(parser.data, data_dir)

    with open(shell_name,"w") as f:
        ## CNVTasks.PreprocessIntervals
        f.write("""{gatk} --java-options \"-Xms{java_mem}G\" \\
            PreprocessIntervals \\
            -L {intervals} \\
            --reference {ref_fasta} \\
            --padding 250 \\
            --bin-length 0 \\
            --interval-merging-rule OVERLAPPING_ONLY \\
            --output targets.preprocessed.interval_list \n""".format(**config_dict))
        # CNVTasks.CollectCounts as CollectCountsTumor
        # CNVTasks.CollectCounts as CollectCountsNormal
        f.write("""{gatk} --java-options \"-Xmx{java_mem}G\"  \\
            CollectReadCounts \\
            -L targets.preprocessed.interval_list \\
            --reference {ref_fasta} \\
            --input  {tumor_bam} \\
            --format HDF5 \\
            --interval-merging-rule OVERLAPPING_ONLY \\
            --output {sample_name}.Tumor.counts.hdf5 \n""".format(**config_dict))
        f.write("""{gatk} --java-options \"-Xmx{java_mem}G\"  \\
            CollectReadCounts \\
            -L targets.preprocessed.interval_list \\
            --reference {ref_fasta} \\
            --input  {normal_bam} \\
            --format HDF5 \\
            --interval-merging-rule OVERLAPPING_ONLY \\
            --output {sample_name}.Normal.counts.hdf5 \n""".format(**config_dict))
        # CNVTasks.CollectAllelicCounts as CollectAllelicCountsTumor
        # CNVTasks.CollectAllelicCounts as CollectAllelicCountsNormal
        f.write("""{gatk} --java-options \"-Xmx{java_mem}G\"  \\
            CollectAllelicCounts \\
            -I {tumor_bam}  \\
            -R {ref_fasta} \\
            -L {common_sites} \\
            --minimum-base-quality 20 \\
            -O {sample_name}.Tumor.allelicCounts.tsv \n""".format(**config_dict))
        f.write("""{gatk} --java-options \"-Xmx{java_mem}G\"  \\
            CollectAllelicCounts \\
            -I {normal_bam} \\
            -R {ref_fasta} \\
            -L {common_sites} \\
            --minimum-base-quality 20 \\
            -O {sample_name}.Normal.allelicCounts.tsv \n""".format(**config_dict))
        # DenoiseReadCounts as DenoiseReadCountsTumor
        # DenoiseReadCounts as DenoiseReadCountsNormal
        f.write("""{gatk} --java-options \"-Xmx{java_mem}G\"  \\
            DenoiseReadCounts \\
            -I {sample_name}.Tumor.readcounts  \\
            --count-panel-of-normals ~{read_count_pon} \\
            --standardized-copy-ratios ~{sample_name}.Tumor.standardizedCR.tsv \\
            --minimum-base-quality 20 \\
            --denoised-copy-ratios ~{sample_name}.Tumor.denoisedCR.tsv \n""".format(**config_dict))
        f.write("""{gatk} --java-options \"-Xmx{java_mem}G\"  \\
            DenoiseReadCounts \\
            -I {sample_name}.Normal.readcounts  \\
            --count-panel-of-normals ~{read_count_pon} \\
            --standardized-copy-ratios ~{sample_name}.Normal.standardizedCR.tsv \\
            --minimum-base-quality 20 \\
            --denoised-copy-ratios ~{sample_name}.Normal.denoisedCR.tsv \n""".format(**config_dict))
        # ModelSegments as ModelSegmentsTumor
        # ModelSegments as ModelSegmentsNormal
        f.write("""{gatk} --java-options \"-Xmx{java_mem}G\"  \\
            ModelSegments \\
            --denoised-copy-ratios ~Tumor.denoised_copy_ratios \\
            --allelic-counts ~Tumor.allelic_counts \\
            --normal-allelic-counts  ollectAllelicCountsNormal.allelic_counts \\
            --minimum-total-allele-count-case ~min_total_allele_count_ \\
            --window-size 0 \\
            --output ~output_dir_Tumor \\
            --output-prefix ~Tumor.entity_id \n""".format(**config_dict))
        f.write("""{gatk} --java-options \"-Xmx{java_mem}G\"  \\
            ModelSegments \\
            --denoised-copy-ratios ~Normal.denoised_copy_ratios \\
            --allelic-counts ~Normal.allelic_counts \\
            --minimum-total-allele-count-case ~min_total_allele_count_ \\
            --window-size 0 \\
            --output ~output_dir_Tumor \\
            --output-prefix ~Normal.entity_id \n""".format(**config_dict))
        
        # CallCopyRatioSegments as CallCopyRatioSegmentsTumor
        # CallCopyRatioSegments as CallCopyRatioSegmentsNormal
        f.write("""{gatk} --java-options \"-Xmx{java_mem}G\"  \\
            CallCopyRatioSegments \\
            --input ~Tumor.copy_ratio_segments \\
            --output ~entity_id.called.seg \n""".format(**config_dict))
        f.write("""{gatk} --java-options \"-Xmx{java_mem}G\"  \\
            CallCopyRatioSegments \\
            --input ~Normal.copy_ratio_segments \\
            --output ~Normal.entity_id.called.seg \n""".format(**config_dict))
        
        # PlotDenoisedCopyRatios as PlotDenoisedCopyRatiosTumor
        # PlotDenoisedCopyRatios as PlotDenoisedCopyRatiosNormal
        f.write("""{gatk} --java-options \"-Xmx{java_mem}G\"  \\
            PlotDenoisedCopyRatios \\
            --standardized-copy-ratios ~DenoiseReadCountsTumor.standardized_copy_ratios \\
            --denoised-copy-ratios ~DenoiseReadCountsTumor.denoised_copy_ratios \\
            --sequence-dictionary {ref_fasta}.dict \\
            --output ./ \\
            --output-prefix CollectCountsTumor.entity_id \n""".format(**config_dict))
        f.write("""{gatk} --java-options \"-Xmx{java_mem}G\"  \\
            PlotDenoisedCopyRatios \\
            --standardized-copy-ratios DenoiseReadCountsNormal.standardized_copy_ratios \\
            --denoised-copy-ratios DenoiseReadCountsNormal.denoised_copy_ratios \\
            --sequence-dictionary {ref_fasta}.dict \\
            --output ./ \\
            --output-prefix ~CollectCountsNormal.entity_id \n""".format(**config_dict))
        
        # PlotModeledSegments as PlotModeledSegmentsTumor
        # PlotModeledSegments as PlotModeledSegmentsNormal
        f.write("""{gatk} --java-options \"-Xmx~{java_mem}G\" \\
            PlotModeledSegments \\
            --denoised-copy-ratios DenoiseReadCountsTumor.denoised_copy_ratios \\
            --allelic-counts ModelSegmentsTumor.het_allelic_counts \\
            --segments ModelSegmentsTumor.modeled_segments \\
            --sequence-dictionary {ref_fasta}.dict \\
            --output ./ \\
            --output-prefix ~CollectCountsTumor.entity_id \n""".format(**config_dict))
        f.write("""{gatk} --java-options \"-Xmx{java_mem}G\" \\
            PlotModeledSegments \\
            --denoised-copy-ratios ~DenoiseReadCountsNormal.denoised_copy_ratios \\
            --allelic-counts ~ModelSegmentsNormal.het_allelic_counts \\
            --segments ~ModelSegmentsNormal.modeled_segments \\
            --sequence-dictionary {ref_fasta}.dict \\
            --output ./ \\
            --output-prefix CollectCountsNormal.entity_id \n""".format(**config_dict))

        # # optional
        # CNVOncotator.CNVOncotatorWorkflow as CNVOncotatorWorkflow

        # CNVFuncotateSegments.CNVFuncotateSegmentsWorkflow as CNVFuncotateSegmentsWorkflow
        # https://github.com/broadinstitute/gatk/blob/master/scripts/cnv_wdl/somatic/cnv_somatic_funcotate_seg_workflow.wdl
        '''
        f.write(""""gatk --java-options "-Xmx~{command_mem_mb}m" FuncotateSegments \\
             --data-sources-path $DATA_SOURCES_FOLDER \\
             --ref-version ~{funcotator_ref_version} \\
             --output-file-format SEG \\
             -R ~{ref_fasta} \\
             --segments ~{input_seg_file} \\
             -O ~{basename_input_seg_file}.funcotated.tsv \\
             ~{interval_list_arg} ~{default="" interval_list} \\
             ~{"--transcript-selection-mode " + transcript_selection_mode} \\
             ~{transcript_selection_arg}~{default="" sep=" --transcript-list " transcript_selection_list} \\
             ~{annotation_def_arg}~{default="" sep=" --annotation-default " annotation_defaults} \\
             ~{annotation_over_arg}~{default="" sep=" --annotation-override " annotation_overrides} \\
             ~{excluded_fields_args}~{default="" sep=" --exclude-field " funcotator_excluded_fields} \\
             ~{extra_args_arg} \n""".format(**config_dict))
        '''
    print("all finished!")

if __name__ == '__main__':
    main()



