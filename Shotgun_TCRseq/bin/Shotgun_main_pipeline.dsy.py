#!/usr/bin/python3.8
import os
import sys
import argparse
import configparser


def _argparse():
    parser = argparse.ArgumentParser(description="This is description")
    parser.add_argument('-c', '--config', action='store', dest='config', default='config.txt', help="this is config file")
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
    if not cf.has_section('Config'):
        os.exit("Error: your config file is not correct.")
    
    # read config file
    config_dict = {
        # 'reference database'
        'ref_fasta': cf.get('Config', 'ref_fasta'),
        'gtf_file': cf.get('Config', 'gtf_file'),

        # project name
        'project_name' : cf.get("Config", "project_name"),
        'sample_name'  : cf.get("Config", "sample_name"),
        'RNA_fastqs' : cf.get("Config", "RNA_fastqs"),
        'scripts_dir': cf.get('Config', 'scripts_dir'),
        # alignment
        'STAR': cf.get('Config', "STAR"),
        'STAR_index': cf.get('Config', "STAR_index"),
        'samtools' : cf.get('Config', "samtools"),

        # expression
        'featureCounts': cf.get('Config', "featureCounts"),
    }
    #make directories:
    project_dir = os.path.abspath(".") + '/' + config_dict['project_name']
    make_dir(project_dir)
    print("# Create work directory")

    # generate shell
    shell_name = project_dir + '/work.' + config_dict['project_name'] + '.RNASeq.sh'
    # only open a file so use try:finally to close.

    with open(shell_name,"w") as f:
        # align
        f.write("{Mixcr} -Xmx64G -Xms64G analyze shotgun --align -OsaveOriginalReads=true --species hs --starting-material rna --receptor-type tcr  -r {sample_name}_TRA.report  {tra_fastqs} {sample_name}_TRA.mixcr.out \n".format(**config_dict))
        f.write("{Mixcr} -Xmx64G -Xms64G analyze shotgun --align -OsaveOriginalReads=true --species hs --starting-material rna --receptor-type tcr  -r {sample_name}_TRB.report  {tra_fastqs} {sample_name}_TRB.mixcr.out \n".format(**config_dict))
        f.write("mkdir -p {sample_name}_TRA {sample_name}_TRB \n".format(**config_dict))
        f.write("mixcr exportReadsForClones -s {sample_name}_TRA.mixcr.out.clna {sample_name}_TRA/{sample_name}_TRA.\n".format(**config_dict))
        f.write("mixcr exportReadsForClones -s {sample_name}_TRB.mixcr.out.clna {sample_name}_TRB/{sample_name}_TRB.\n".format(**config_dict))

        f.write("python3 {scripts_dir}/0.read2barcode.py {TRA_Barcode_file} ../data/{sample_name}_S1_L001_R1_001.fastq.gz {sample_name}_TRA\n".format(**config_dict))
        f.write("python3 {scripts_dir}/0.read2barcode.py {TRB_Barcode_file} ../data/{sample_name}_S1_L001_R1_001.fastq.gz {sample_name}_TRB\n".format(**config_dict))
        
        f.write("python3 {scripts_dir}/1.mixcr2barcode.py G125E2L1_TRA.barcode.read.txt G125E2L1_TRA.mixcr.out.clonotypes.TRA.txt TRA\n".format(**config_dict))

        # add TPM and add RNAseq read depth.
        f.write("python3 {scripts_dir}/add_TPM.py {sample_name}.variants.funcotated.with.minigene.MAF.xls {sample_name}.RNAseq.transcript.counts.TPM.txt transcript_id\n".format(**config_dict))
        f.write("python3 {scripts_dir}/add_RNAseq_read_depth.py {sample_name}.variants.funcotated.with.minigene.MAF.add.TPM.xls {sample_name}.RNAseq.mpileup.txt\n".format(**config_dict))
    print("all finished!")


if __name__ == '__main__':
    main()
