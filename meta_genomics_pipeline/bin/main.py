#!/usr/bin/python3
import os
import sys
import argparse
import configparser


def _argparse():
    parser = argparse.ArgumentParser(description="This is a metagenomics pipeline")
    parser.add_argument('-c', '--config', action='store', dest='config', default='config.txt', help = "this is config file.")
    parser.add_argument('-l', '--list', action='store', dest='data_file', default='raw.data.list.xls', help = "your data list file.")
    return parser.parse_args()

def make_dir(*dir):
    for each in dir:
        if not os.path.exists(each):
            os.mkdir(each)

def read_raw_data(data_file):
    data_dict = {}
    with open(data_file, "r") as f:
        for line in f:
            line = line.rstrip("\n")
            sample_name, fastq = line.split()
            #fastq1,fastq2 = fastq.split((",",";"))
            data_dict[sample_name] = fastq
    return data_dict

def main():
    """ docstring for __main__"""
    parser = _argparse()
    cf = configparser.ConfigParser(allow_no_value=True)
    cf.read(parser.config)
    if not cf.has_section('config'):
        sys.exit("Error: your config file is not correct.")

    # read config file:
    config_dict = {
        # script dir
        'script_dir' : cf.get('config', 'script_dir'),
        # reference database
        'kneaddata_ref': cf.get('config', 'kneaddata_ref'),
        'kraken2_db': cf.get('config', 'kraken2_db'),
        # software.
        'kneaddata' : cf.get('config', 'kneaddata'),
        'kraken2' : cf.get('config', 'kraken2'),
        'bracken' : cf.get('config', 'bracken'),
        'kreport2mpa': cf.get('config', 'kreport2mpa'),
        'humann2' : cf.get('config', 'humann2'),
        'megahit' : cf.get('config', 'megahit'),
        'quast' : cf.get('config', 'quast'),
        'prokka' : cf.get('config', 'prokka'),
        'cd-hit' : cf.get('config','cd-hit'),
        'salmon' : cf.get('config','salmon'),
        # project name.
        'project_name' : cf.get('config','project_name'),
        }
    
    project_dir = os.path.abspath(".") + '/' + config_dict['project_name']
    make_dir(project_dir)
    kneaddata_dir = os.path.abspath(project_dir) + '/kneaddata_result'
    taxonomic_dir = os.path.abspath(project_dir) + '/kraken2_result'
    humann2_dir = os.path.abspath(project_dir) +'/humann2_result'
    megahit_dir   = os.path.abspath(project_dir) + '/megahit_result'
    
    config_dict.update({"project_dir":project_dir,"kneaddata_dir":kneaddata_dir,"taxonomic_dir":taxonomic_dir,"humann2_dir":humann2_dir,"megahit_dir":megahit_dir})
    make_dir(kneaddata_dir, taxonomic_dir, humann2_dir, megahit_dir)
    print("# Create work directory")
    # print(config_dict)

    # generate shell
    shell_name = project_dir + '/work.' + config_dict['project_name'] + '.metagenomics.sh'
    # deal fastq:
    fastq_data_dict = read_raw_data(parser.data_file)
    with open(shell_name, "w") as f:
        # step1: filtering
        for each_sample in fastq_data_dict:
            f.write("{kneaddata} --input {0} --reference-db {kneaddata_ref} --output {kneaddata_dir}  --output-prefix {1}.kneaddata \n".format(fastq_data_dict[each_sample], each_sample, **config_dict))
        # step2:
            f.write("{kraken2} --db {kraken2_db}  --threads 50  --report {taxonomic_dir}/{0}.kraken2.report --output {taxonomic_dir}/{0}.kraken2.output {kneaddata_dir}/{0}.kneaddata.fastq 2>{taxonomic_dir}/{0}.kraken2.log\n".format(each_sample, **config_dict))
            f.write("{bracken} -d {kraken2_db} -i {taxonomic_dir}/{0}.kraken2.report -o {taxonomic_dir}/{0}.S.bracken -w {taxonomic_dir}/{0}.S.bracken.report -r 150 -l S \n".format(each_sample, **config_dict))
            f.write("{kreport2mpa} -r {taxonomic_dir}/{0}.S.bracken.report -o {taxonomic_dir}/{0}.S.kreport2mpa.report\n".format(each_sample, **config_dict))

            # step3: see: https://github.com/biobakery/biobakery/wiki/humann2  3.1, 3.2
            f.write("{humann2} --verbose --remove-temp-output --nucleotide-database /cygene/software/biosoftware/metagenomics/humann2_data/chocophlan --threads 50 --input {kneaddata_dir}/{0}.kneaddata.fastq --output {humann2_dir} \n".format(each_sample, **config_dict))
            f.write("{script_dir}/humann2_rename_table --input {humann2_dir}/{0}_genefamilies.tsv --output {humann2_dir}/{0}_genefamilies-names.tsv --names uniref90\n".format(each_sample, **config_dict))
            # Normalizing RPKs to relative abundance
            f.write("{script_dir}/humann2_renorm_table --input {humann2_dir}/{0}_genefamilies.tsv --output {humann2_dir}/{0}_genefamilies-cpm.tsv --units cpm --update-snames\n".format(each_sample, **config_dict))
            # will generate a megahit_dir.
            f.write("{megahit} -r {kneaddata_dir}/{0}.kneaddata.fastq -o {megahit_dir}/{0} --out-prefix {0}.megahit.final \n".format(each_sample, **config_dict))
            f.write("{quast} {megahit_dir}/{0}/{0}.megahit.final.contigs.fa -o {megahit_dir}/{0}/{0}.megahit-quast-report \n".format(each_sample, **config_dict))

            f.write("# \n".format(**config_dict))

            f.write("{prokka} {megahit_dir}/{0}/{0}.megahit.final.contigs.fa --outdir {project_dir}/prokka_annotation --force --prefix {0} --metagenome --kingdom Bacteria \n".format(each_sample, **config_dict))

            f.write("# salmon\n".format(**config_dict))
        f.write("{script_dir}/kneaddata_read_count_table --input {kneaddata_dir} --output Summary.kneadata.results.xls \n".format(**config_dict))
        f.write("{script_dir}/summary_kraken_count_table.py --input {taxonomic_dir} --output Summary.kraken2.results.xls \n".format(**config_dict))
        f.write("{script_dir}/merge_metaphlan_tables.py {humann2_dir}/*/*.tsv > Summary.metaphlan.results.xls\n".format(**config_dict))
        f.write("{script_dir}/humann2_join_tables --input {humann2_dir} --output Summary.humann2.results.xls\n".format(**config_dict))

    print("all finished!")

if __name__ == '__main__':
    main()
