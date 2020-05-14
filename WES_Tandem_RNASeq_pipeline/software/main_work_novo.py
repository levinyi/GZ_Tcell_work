import os
import sys
import argparse
import ConfigParser


def _argparse():
    parser = argparse.ArgumentParser(description="This is description")
    parser.add_argument('-c', '--config', action='store', dest='config', default='config.txt', help="this is config file")
    parser.add_argument('-l', '--list', action='store', dest='data', default='data.xls', help="this is data list file")
    return parser.parse_args()


def make_dir(*dir):
    for each in dir:
        if not os.path.exists(each):
            os.mkdir(each, 0755)



def main():
    """docstring for __main__"""
    parser = _argparse()
    cf = ConfigParser.ConfigParser(allow_no_value=True)
    cf.read(parser.config)
    excel_file = parser.data
    excel_file = os.path.abspath(excel_file)
    # check config file
    if not cf.has_section('Config'):
        os.exit("Error: your config file is not correct.")
    
    # read config file
    config_dict = {
        'python' : cf.get('Config', 'python'),
        'software_path' : cf.get('Config', 'software_path'),
        'database' : cf.get('Config', 'database'),
        'project_name' : cf.get("Config", "project_name"),
        'RNAfastqR1' : cf.get('Config',"RNASeqFastqData").split(",")[0],
        'RNAfastqR2' : cf.get('Config',"RNASeqFastqData").split(",")[1],
        'WESfastqR1' : cf.get('Config',"WESSeqFastqData").split(",")[0],
        'WESfastqR2' : cf.get('Config',"WESSeqFastqData").split(",")[1],

        #'WES_bam_data' : cf.get("Config", "WESBAMData"),
        'excel_file' : excel_file,
    }
    #make directories:
    project_dir = os.path.abspath(".") + '/' + config_dict['project_name']
    make_dir(project_dir)
    print("# Create work directory")

    # generate shell
    shell_name = project_dir + '/work.' + config_dict['project_name'] + '.sh'
    #shell_name = shell_dir + '/work.' + config_dict['project_name'] + '.sh'
    # only open a file so use try:finally to close.
    # rawdata_dict = deal_rawdata(parser.data, data_dir)

    with open(shell_name,"w") as f:
        # deal with the standard input File format.
        f.write("#{python} {software_path}/deal_input.STF.py {excel_file} > {project_name}.STF \n".format(**config_dict))
        # replace mutation position in genome file.
        f.write("#{python} {software_path}/replace_genome.novo.py {database}/hg19.fa {project_name}.STF indel.info.txt > genome_replaced.info.xls \n".format(**config_dict))
        f.write("#echo \"finished step:genome_replace\"\n\n")
        f.write("#if [ -f genome_replaced_mutation.fa.fai ]; then\n")
        f.write("#\t rm genome_replaced_mutation.fa.fai\n".format(**config_dict))
        f.write("#fi\n\n")
        f.write("#{software_path}/gffread {database}/hg19.filtered.gtf -g genome_replaced_mutation.fa -x genome_replaced_mutation.cds.fa\n\n".format(**config_dict))
        f.write("#echo \"finished step: gffread\"\n\n")

        # check the replaced fa.
        f.write("#{python} {software_path}/compare_Two_Fasta.py --ref {database}/hg19.fa --treat genome_replaced_mutation.fa --flag long -o genome.VS.replaced.comparing.txt\n".format(**config_dict))
        f.write("#{python} {software_path}/compare_Two_Fasta.py --ref {database}/hg19.cds.fa --treat genome_replaced_mutation.cds.fa --flag short -o genome_CDS.VS.replaced.comparing.txt\n".format(**config_dict))
        f.write("#{python} {software_path}/check_replaced.report.py --ref genome.VS.replaced.comparing.txt --stf {project_name}.STF --flag genome -o genome_mutation.checked.out\n".format(**config_dict))
        f.write("#{python} {software_path}/check_replaced.report.py --ref genome_CDS.VS.replaced.comparing.txt --stf {project_name}.STF --flag CDS -o genome_CDS_mutation.checked.out\n".format(**config_dict))

        f.write("#if [[ ! -s  indel.info.txt ]] ; then \n")  # no insert info.
        # remove N in cds
        f.write("#\t {python} {software_path}/deal.CDS.insert.py -c genome_replaced_mutation.cds.fa --gtf {database}/hg19.filtered.gtf -o genome_replaced_mutation.rmN.cds.fa\n".format(**config_dict))
        f.write("#else \n") # insert info.
        f.write("#\t {python} {software_path}/deal.CDS.insert.py -c genome_replaced_mutation.cds.fa --gtf {database}/hg19.filtered.gtf -l indel.info.txt -o genome_replaced_mutation.rmN.cds.fa\n".format(**config_dict))
        f.write("#\t echo 'finished deal.cds.insert.py' \n")
        f.write("#fi\n\n")

        f.write("#{software_path}/transeq genome_replaced_mutation.rmN.cds.fa genome_replaced_mutation.rmN.pep.fa\n\n".format(**config_dict))
        f.write("{python} {software_path}/compare_Two_Fasta.py --ref {database}/hg19.pep.fa --treat genome_replaced_mutation.rmN.pep.fa -o genome_pep.VS.replaced.comparing.txt\n".format(**config_dict))
        f.write("#{python} {software_path}/check_replaced.report.py --ref genome_pep.VS.replaced.comparing.txt --stf {project_name}.STF --flag PEP -o genome_PEP_mutation.checked.out\n".format(**config_dict))

        f.write("#{python} {software_path}/check.pep.py {database}/hg19.pep.fa {project_name}.STF genome_replaced_mutation.rmN.pep.fa >{project_name}.firstStep.Mutated.Tandem.Minigenes.xls\n".format(**config_dict))
        f.write("#echo \"finished step: check pep\"\n\n")

        f.write("{python} {software_path}/extract_bed.py {project_name}.firstStep.Mutated.Tandem.Minigenes.xls >{project_name}.snp.checked.bed\n".format(**config_dict))
        f.write("echo \"finished step:extract bed\"\n\n")
        # step 1 deal RNAseq data.
        f.write("# {software_path}/STAR --runThreadN 10 --genomeDir {database}/hg19_star2.7_index  --readFilesCommand zcat --readFilesIn {RNAfastqR1}  {RNAfastqR2}  --outFileNamePrefix {project_name}.RNA.  --outSAMtype BAM SortedByCoordinate \n".format(**config_dict ))
        f.write("# echo \"finished step: STAR RNA\"\n\n")
        f.write("# {software_path}/samtools index {project_name}.RNA.Aligned.sortedByCoord.out.bam \n".format(**config_dict))
        f.write("# echo \"finished step: samtools index\"\n\n")
        f.write("{software_path}/samtools mpileup -l {project_name}.snp.checked.bed -f {database}/hg19.fa {project_name}.RNA.Aligned.sortedByCoord.out.bam  >{project_name}.RNAseq.mpileup.txt \n".format(**config_dict))
        f.write("echo \"finished step: samtools mpileup\"\n\n")
        
        f.write("# {software_path}/featureCounts -O -T 20 -t exon -g gene_id -a {database}/hg19.filtered.gtf -o {project_name}.gene.counts.txt  {project_name}.RNA.Aligned.sortedByCoord.out.bam\n".format(**config_dict))
        f.write("# {python} {software_path}/featureCounts2TPM.py -a {project_name}.gene.counts.txt -o {project_name}.RNAseq.gene.counts.TPM.txt\n".format(**config_dict))
        f.write("# echo \"finished step: calculate TPM\"\n\n")

        # step 2 WES dat
        f.write("# {software_path}/STAR --runThreadN 10 --genomeDir {database}/hg19_star2.7_index  --readFilesCommand zcat --readFilesIn {WESfastqR1}  {WESfastqR2}  --outFileNamePrefix {project_name}.WES.  --outSAMtype BAM SortedByCoordinate \n".format(**config_dict ))
        f.write("# {software_path}/samtools index {project_name}.WES.Aligned.sortedByCoord.out.bam \n".format(**config_dict))
        f.write("{software_path}/samtools mpileup -l {project_name}.snp.checked.bed -f {database}/hg19.fa {project_name}.WES.Aligned.sortedByCoord.out.bam  >{project_name}.WES.mpileup.txt \n".format(**config_dict))
        f.write("echo \"finished step: samtools mpileup WES\"\n\n")

        # step3 combin results:
        f.write("{python} {software_path}/combine_WES_RNA_result.py {project_name}.firstStep.Mutated.Tandem.Minigenes.xls {project_name}.RNAseq.gene.counts.TPM.txt {project_name}.RNAseq.mpileup.txt {project_name}.WES.mpileup.txt > {project_name}.TableS4.Mutated.Tandem.Minigenes.V4.xls \n".format(**config_dict))
        f.write("echo \"finished step: combine\"\n\n")
    print("all finished!")

if __name__ == '__main__':
    main()



