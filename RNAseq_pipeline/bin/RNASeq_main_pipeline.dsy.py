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
        f.write("{STAR} --runThreadN 20 --genomeDir {STAR_index} --readFilesCommand zcat --readFilesIn {RNA_fastqs}  --outFileNamePrefix {sample_name}.RNA. --quantMode GeneCounts --outSAMtype BAM SortedByCoordinate \n".format(**config_dict))
        f.write("{samtools} index {sample_name}.RNA.Aligned.sortedByCoord.out.bam\n".format(**config_dict))
        f.write("{samtools depth -@ 20 -b {sample_name}.snp.checked.bed -a {sample_name}.RNA.Aligned.sortedByCoord.out.bam -o {sample_name}.RNAseq.snp.checked.depth.txt\n}".format(**config_dict))
        f.write("{samtools} mpileup -ABQ0 -l {sample_name}.snp.checked.bed -f {ref_fasta} {sample_name}.RNA.Aligned.sortedByCoord.out.bam -aa -o {sample_name}.RNAseq.mpileup.txt\n".format(**config_dict))

        f.write("{featureCounts} -O -T 20 -t exon -g gene_name -p -a {gtf_file} -o {sample_name}.RNAseq.exon.gene_name.counts.txt  {sample_name}.RNA.Aligned.sortedByCoord.out.bam \n".format(**config_dict))
        f.write("python3 {scripts_dir}/featureCounts2TPM.py -a {sample_name}.RNAseq.exon.gene_name.counts.txt -o {sample_name}.RNAseq.exon.gene_name.counts.TPM.txt\n".format(**config_dict))

        # add TPM and add RNAseq read depth.
        f.write("python3 {scripts_dir}/add_TPM.py {sample_name}.variants.funcotated.with.minigene.MAF.xls {sample_name}.RNAseq.exon.gene_name.counts.TPM.txt gene_id\n".format(**config_dict))
        f.write("python3 {scripts_dir}/add_RNAseq_read_depth.py {sample_name}.variants.funcotated.with.minigene.MAF.add.TPM.xls {sample_name}.RNAseq.mpileup.txt\n".format(**config_dict))
    print("all finished!")


if __name__ == '__main__':
    main()
