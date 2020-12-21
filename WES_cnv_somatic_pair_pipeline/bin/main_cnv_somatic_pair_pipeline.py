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
        'funcotator_ref_version' : cf.get("GATK", "funcotator_ref_version"),
        
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
            --output targets.preprocessed.interval_list \n\n""".format(**config_dict))

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
            --output {sample_name}.Normal.counts.hdf5 \n\n""".format(**config_dict))

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
        f.write("""{gatk} --java-options \"-Xmx{java_mem}G\"  \\
        CreateReadCountPanelOfNormals \\
            -I {sample_name}.Normal.counts.hdf5  \\
            -O cnv.pon.hdf5 \n""".format(**config_dict))
        # DenoiseReadCounts as DenoiseReadCountsTumor
        # DenoiseReadCounts as DenoiseReadCountsNormal
        f.write("""{gatk} --java-options \"-Xmx{java_mem}G\"  \\
        DenoiseReadCounts \\
            -I {sample_name}.Tumor.counts.hdf5  \\
            --count-panel-of-normals cnv.pon.hdf5 \\
            --standardized-copy-ratios {sample_name}.Tumor.standardizedCR.tsv \\
            --denoised-copy-ratios {sample_name}.Tumor.denoisedCR.tsv \n""".format(**config_dict))
        f.write("""{gatk} --java-options \"-Xmx{java_mem}G\"  \\
        DenoiseReadCounts \\
            -I {sample_name}.Normal.counts.hdf5  \\
            --count-panel-of-normals cnv.pon.hdf5 \\
            --standardized-copy-ratios {sample_name}.Normal.standardizedCR.tsv \\
            --denoised-copy-ratios {sample_name}.Normal.denoisedCR.tsv \n\n""".format(**config_dict))

        # ModelSegments as ModelSegmentsTumor
        # ModelSegments as ModelSegmentsNormal
        f.write("""{gatk} --java-options \"-Xmx{java_mem}G\"  \\
        ModelSegments \\
            --denoised-copy-ratios {sample_name}.Tumor.denoisedCR.tsv \\
            --allelic-counts {sample_name}.Tumor.allelicCounts.tsv \\
            --normal-allelic-counts {sample_name}.Normal.allelicCounts.tsv \\
            --output ./ \\
            --output-prefix {sample_name}.Tumor \n""".format(**config_dict))
        f.write("""{gatk} --java-options \"-Xmx{java_mem}G\"  \\
        ModelSegments \\
            --denoised-copy-ratios {sample_name}.Normal.denoisedCR.tsv \\
            --allelic-counts {sample_name}.Normal.allelicCounts.tsv \\
            --output ./ \\
            --output-prefix {sample_name}.Normal \n\n""".format(**config_dict))
        
        # CallCopyRatioSegments as CallCopyRatioSegmentsTumor
        # CallCopyRatioSegments as CallCopyRatioSegmentsNormal
        f.write("""{gatk} --java-options \"-Xmx{java_mem}G\"  \\
        CallCopyRatioSegments \\
            --input {sample_name}.Tumor.cr.seg \\
            --output {sample_name}.Tumor.called.seg \n""".format(**config_dict))
        f.write("""{gatk} --java-options \"-Xmx{java_mem}G\"  \\
        CallCopyRatioSegments \\
            --input {sample_name}.Normal.cr.seg \\
            --output {sample_name}.Normal.called.seg \n\n""".format(**config_dict))
        
        # PlotDenoisedCopyRatios as PlotDenoisedCopyRatiosTumor
        # PlotDenoisedCopyRatios as PlotDenoisedCopyRatiosNormal
        f.write("""{gatk} --java-options \"-Xmx{java_mem}G\"  \\
        PlotDenoisedCopyRatios \\
            --standardized-copy-ratios {sample_name}.Tumor.standardizedCR.tsv  \\
            --denoised-copy-ratios {sample_name}.Tumor.denoisedCR.tsv  \\
            --sequence-dictionary {ref_fasta}.dict \\
            --output ./ \\
            --output-prefix {sample_name}.Tumor \n""".format(**config_dict))
        f.write("""{gatk} --java-options \"-Xmx{java_mem}G\"  \\
        PlotDenoisedCopyRatios \\
            --standardized-copy-ratios {sample_name}.Normal.standardizedCR.tsv \\
            --denoised-copy-ratios {sample_name}.Normal.denoisedCR.tsv \\
            --sequence-dictionary {ref_fasta}.dict \\
            --output ./ \\
            --output-prefix {sample_name}.Normal \n\n""".format(**config_dict))

        # PlotModeledSegments as PlotModeledSegmentsTumor
        # PlotModeledSegments as PlotModeledSegmentsNormal
        f.write("""{gatk} \\
        PlotModeledSegments \\
            --denoised-copy-ratios {sample_name}.Tumor.denoisedCR.tsv \\
            --allelic-counts {sample_name}.Tumor.hets.tsv \\
            --segments {sample_name}.Tumor.modelFinal.seg \\
            --sequence-dictionary {ref_fasta}.dict \\
            --output ./ \\
            --output-prefix {sample_name}.Tumor \n""".format(**config_dict))
        f.write("""{gatk} --java-options \"-Xmx{java_mem}G\" \\
        PlotModeledSegments \\
            --denoised-copy-ratios {sample_name}.Normal.denoisedCR.tsv \\
            --allelic-counts {sample_name}.Normal.hets.tsv \\
            --segments {sample_name}.Normal.modelFinal.seg \\
            --sequence-dictionary {ref_fasta}.dict \\
            --output ./ \\
            --output-prefix {sample_name}.Normal \n\n""".format(**config_dict))

        # CNVFuncotateSegments.CNVFuncotateSegmentsWorkflow as CNVFuncotateSegmentsWorkflow
        # https://github.com/broadinstitute/gatk/blob/master/scripts/cnv_wdl/somatic/cnv_somatic_funcotate_seg_workflow.wdl
        f.write("""{gatk} --java-options \"-Xmx{java_mem}G\"  \\
        FuncotateSegments \\
            --data-sources-path {funcotator_dataSources} \\
            --ref-version {funcotator_ref_version} \\
            --output-file-format SEG \\
            -R {ref_fasta} \\
            --segments {sample_name}.Normal.modelFinal.seg \\
            -O {sample_name}.Normal.modelFinal.seg.funcotated.tsv \n\n""".format(**config_dict))
        f.write("""{gatk} --java-options \"-Xmx{java_mem}G\"  \\
        FuncotateSegments \\
            --data-sources-path {funcotator_dataSources} \\
            --ref-version {funcotator_ref_version} \\
            --output-file-format SEG \\
            -R {ref_fasta} \\
            --segments {sample_name}.Tumor.modelFinal.seg \\
            -O {sample_name}.Tumor.modelFinal.seg.funcotated.tsv \n\n""".format(**config_dict))
    print("all finished!")

if __name__ == '__main__':
    main()
