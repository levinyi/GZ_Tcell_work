#!/usr/bin/python3.8
import os
import argparse
import configparser
import logging


def _argparse():
    parser = argparse.ArgumentParser(description="This is description")
    parser.add_argument('-c', '--config', action='store', dest='config', default='config.txt', help="this is config file")
    return parser.parse_args()


def make_dir(*dir):
    for each in dir:
        if not os.path.exists(each):
            os.makedirs(each)


def main():
    """docstring for main"""
    logging.basicConfig(level=logging.DEBUG, format="%(asctime)s %(levelname)s %(message)s", datefmt='%Y-%m-%d %H:%M:%S')
    parser = _argparse()

    cf = configparser.ConfigParser(allow_no_value=True)
    cf.read(parser.config)

    # check config file
    if not cf.has_section('Config'):
        logging.error("\033[1;31mError: \033[0m Your config file is not correct.")
    
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
    logging.info("Created working directory")

    # generate shell
    shell_name = project_dir + '/work.' + config_dict['project_name'] + '.RNASeq2TMP.sh'
    # only open a file so use try:finally to close.

    with open(shell_name,"w") as f:
        # align
        f.write("{STAR} --runThreadN 10 --genomeDir {STAR_index} --readFilesCommand zcat --readFilesIn {RNA_fastqs}  --outFileNamePrefix {sample_name}.RNA. --outSAMtype BAM SortedByCoordinate \n".format(**config_dict))
        f.write("{samtools} index {sample_name}.RNA.Aligned.sortedByCoord.out.bam\n".format(**config_dict))

        f.write("{featureCounts} -O -T 20 -t exon -g gene_name -p -a {gtf_file} -o {sample_name}.RNAseq.exon.gene_name.counts.txt  {sample_name}.RNA.Aligned.sortedByCoord.out.bam \n".format(**config_dict))
        f.write("python3 {scripts_dir}/featureCounts2TPM.py -a {sample_name}.RNAseq.exon.gene_name.counts.txt -o {sample_name}.RNAseq.exon.gene_name.counts.TPM.txt\n".format(**config_dict))
    logging.info("\033[0;36mSuccessful!\033[0m")


if __name__ == '__main__':
    main()
