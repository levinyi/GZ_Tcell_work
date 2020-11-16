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
        # reference database
        'kneaddata_ref': cf.get('config', 'kneaddata_ref'),
        'kraken2_db': cf.get('config', 'kraken2_db')
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
    taxonomic_dir = os.path.abspath(project_dir) + '/kraken2'
    functional_dir = os.path.abspath(project_dir) + '/functional_result'
    megahit_dir   = os.path.abspath(project_dir) + '/megahit_result'
    
    config_dict.update({"project_dir":project_dir,"kneaddata_dir":kneaddata_dir,"taxonomic_dir":taxonomic_dir,"functional_dir":functional_dir,"megahit_dir":megahit_dir})
    make_dir(kneaddata_dir, taxonomic_dir, functional_dir, megahit_dir)
    print("# Create work directory")
    # print(config_dict)

    # generate shell
    shell_name = project_dir + '/work.' + config_dict['project_name'] + '.metagenomics.sh'
    # deal fastq:
    fastq_data_dict = read_raw_data(parser.data_file)
    with open(shell_name, "w") as f:
        # step1: filtering
        for each_sample in fastq_data_dict:
            f.write("{kneaddata} --input {} --reference-db {kneaddata_ref} --output {kneaddata_dir}\n".format(fastq_data_dict[each_sample],**config_dict))
        # step2:
            f.write("{kraken2} --db {kraken2_db}  --threads 50  --report ./{}.report --output ./{}.output  {}eaddata.fastq 2>log\n".format(**config_dict))
            f.write("{bracken} -d {kraken2_db} -i {}.report -o {}.S.bracken -w {}.S.bracken.report -r 150 -l S \n".format(each_sample, **config_dict))
            f.write("{kreport2mpa} -r {}.S.bracken.report   -o {}.new.report\n".format(each_sample, **config_dict))

            # step3:
            f.write("humann2\n".format(**config_dict))

            f.write("{megahit} -r {}.kneaddata.fastq -o {} --out-prefix {}.megahit.final \n".format(each_sample, **config_dict))
            f.write("{quast} {outdir}/{}.megahit.final.contigs.fa -o combined-report \n".format(**config_dict))

            f.write("metaquast\n".format(**config_dict))

            f.write("prokka\n".format(**config_dict))

            f.write("salmon\n".format(**config_dict))

    print("all finished!")

if __name__ == '__main__':
    main()
