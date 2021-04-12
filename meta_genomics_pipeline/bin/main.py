#!/usr/bin/python3
import os
import sys
import argparse
import configparser


def usage():
    print('''
    updates:
    
    20210411:   add PE data parsing.
    20200611:   created.
    '''.format())

def _argparse():
    parser = argparse.ArgumentParser(description="This is metagenomics pipeline")
    parser.add_argument('-c', '--config', action='store', dest='config', default='config.txt', help = "this is config file.")
    parser.add_argument('-l', '--list', action='store', dest='data_file', default='data.list', help = "your data list file.")
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
        'humann2_db': cf.get('config', 'humann2_db'),
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
    fastqc_dir = os.path.abspath(project_dir) + '/fastqc_result'
    kneaddata_dir = os.path.abspath(project_dir) + '/kneaddata_result'
    taxonomic_dir = os.path.abspath(project_dir) + '/kraken2_result'
    humann2_dir = os.path.abspath(project_dir) +'/humann2_result'
    prokka_dir  = os.path.abspath(project_dir) + '/prokka_result'
    megahit_dir   = os.path.abspath(project_dir) + '/megahit_result'
    summary_dir   = os.path.abspath(project_dir) + '/summary_report'
    salmon_dir    = os.path.abspath(project_dir) + '/salmon_result'
    
    config_dict.update(
            {"project_dir" : project_dir,
              "fastqc_dir" : fastqc_dir,
           "kneaddata_dir" : kneaddata_dir,
           "taxonomic_dir" : taxonomic_dir,
             "humann2_dir" : humann2_dir,
             "megahit_dir" : megahit_dir, 
             "salmon_dir"  : salmon_dir,
             "prokka_dir"  : prokka_dir,
             })
    make_dir(fastqc_dir, kneaddata_dir, taxonomic_dir, humann2_dir, megahit_dir, salmon_dir, summary_dir)
    print("# Create work directory")
    # print(config_dict)

    # generate shell
    shell_name = project_dir + '/work.' + config_dict['project_name'] + '.metagenomics.sh'
    # deal fastq:
    fastq_data_dict = read_raw_data(parser.data_file)
    # print(fastq_data_dict.values())
    ### fastq_data_dict :{'sample01': 'houmaohu_R1.fq.gz', 'sample02': '/YJ2021620_R1.fq.gz'}
    with open(shell_name, "w") as f:
        # step1: filtering.
        for each_sample in fastq_data_dict:
            f.write("ln -s {0} {fastqc_dir}/{1}.raw.fastq.gz \n".format(fastq_data_dict[each_sample], each_sample, **config_dict))
            f.write("echo \"finshed link fastq files\" \n")
            f.write("{kneaddata} --input {0} --reference-db {kneaddata_ref} --output {kneaddata_dir} --output-prefix {1}.kneaddata --run-trim-repetitive --sequencer-source none -t 10 -p 10 --fastqc fastqc \n".format(fastq_data_dict[each_sample], each_sample, **config_dict))
            f.write("echo \"finished kneaddata.\" \n")
        # step2:
            f.write("{kraken2} --db {kraken2_db}  --threads 50  --report {taxonomic_dir}/{0}.kraken2.report --output {taxonomic_dir}/{0}.kraken2.output {kneaddata_dir}/{0}.kneaddata.fastq 2>{taxonomic_dir}/{0}.kraken2.log \n".format(each_sample, **config_dict))
            f.write("{bracken} -d {kraken2_db} -i {taxonomic_dir}/{0}.kraken2.report -o {taxonomic_dir}/{0}.S.bracken -w {taxonomic_dir}/{0}.S.bracken.report -r 150 -l S \n".format(each_sample, **config_dict))
            f.write("{kreport2mpa} -r {taxonomic_dir}/{0}.S.bracken.report -o {taxonomic_dir}/{0}.S.kreport2mpa.report \n".format(each_sample, **config_dict))

            # step3: see: https://github.com/biobakery/biobakery/wiki/humann2  3.1, 3.2
            f.write("{humann2} --verbose --remove-temp-output --nucleotide-database {humann2_db}/chocophlan --threads 50 --input {kneaddata_dir}/{0}.kneaddata.fastq --output {humann2_dir} \n".format(each_sample, **config_dict))
            f.write("{script_dir}/humann2_rename_table --input {humann2_dir}/{0}_genefamilies.tsv --output {humann2_dir}/{0}_genefamilies-names.tsv --names uniref90 \n".format(each_sample, **config_dict))
            # Normalizing RPKs to relative abundance
            f.write("{script_dir}/humann2_renorm_table --input {humann2_dir}/{0}_genefamilies.tsv --output {humann2_dir}/{0}_genefamilies-cpm.tsv --units cpm --update-snames \n".format(each_sample, **config_dict))
            # will generate a megahit_dir.
            f.write("{megahit} -r {kneaddata_dir}/{0}.kneaddata.fastq -o {megahit_dir}/{0} --out-prefix {0}.megahit.final \n".format(each_sample, **config_dict))
            f.write("{quast} {megahit_dir}/{0}/{0}.megahit.final.contigs.fa -o {megahit_dir}/{0}/{0}.megahit-quast-report \n".format(each_sample, **config_dict))

            f.write("# \n".format(**config_dict))

            f.write("{prokka} {megahit_dir}/{0}/{0}.megahit.final.contigs.fa --outdir {prokka_dir} --force --prefix {0} --metagenome --kingdom Bacteria --quiet \n".format(each_sample, **config_dict))

            f.write("salmon index --index {salmon_dir}/{0} --transcripts {prokka_dir}/{0}.ffn \n".format(each_sample, **config_dict))
            f.write("salmon quant --index {salmon_dir}/{0} --libType IU -r {kneaddata_dir}/{0}.kneaddata.fastq -o {salmon_dir}/{0} --quiet \n".format(each_sample, **config_dict))
            f.write("python {script_dir}/gather-counts.py {salmon_dir}\n".format(**config_dict))
            # circos
            f.write("# \n")
            f.write("# circos ./circos.config \n".format(**config_dict))

        # fastqc
        f.write("\n# fastqc raw data and filgered data \n")
        f.write("fastqc {fastqc_dir}/*.fastq.gz --outdir {fastqc_dir} -t 6 --quiet \n".format(**config_dict))
        f.write("fastqc {kneaddata_dir}/*.kneaddata.fastq --outdir {fastqc_dir} -t 6 --quiet \n".format(**config_dict))
        f.write("multiqc {fastqc_dir} --outdir {fastqc_dir} --quiet --no-ansi \n".format(**config_dict))
        f.write("\n# summary statistics: \n")

        f.write("{script_dir}/kneaddata_read_count_table --input {kneaddata_dir} --output ./summary_report/Summary.kneadata.results.xls \n".format(**config_dict))
        f.write("{script_dir}/summary_kraken_count_table.py --input {taxonomic_dir} --output ./summary_report/Summary.kraken2.results.xls \n".format(**config_dict))
        # f.write("{script_dir}/merge_metaphlan_tables.py {humann2_dir}/*/*.tsv > ./summary_report/Summary.metaphlan.results.xls \n".format(**config_dict))
        f.write("{script_dir}/humann2_join_tables --input {humann2_dir} --output ./summary_report/Summary.humann2.results.xls \n".format(**config_dict))
        # Join HUMAnN2 output per sample into one table.
        f.write("{script_dir}/humann2_join_tables -s --input {humann2_dir} --file_name pathabundance --output ./summary_report/Summary.humann2.pathabundance.tsv\n".format(**config_dict))
        f.write("{script_dir}/humann2_join_tables -s --input {humann2_dir} --file_name pathcoverage --output ./summary_report/Summary.humann2.pathcoverage.tsv \n".format(**config_dict))
        f.write("{script_dir}/humann2_join_tables -s --input {humann2_dir} --file_name genefamilies --output ./summary_report/Summary.humann2.genefamilies.tsv \n".format(**config_dict))
        # Re-normalize gene family and pathway abundances (so that all samples are in units of copies per million).
        f.write("{script_dir}/humann2_renorm_table --input ./summary_report/Summary.humann2.pathabundance.tsv --units cpm --output ./summary_report/Summary.humann2.pathabundance_cpm.tsv\n".format(**config_dict))
        f.write("{script_dir}/humann2_renorm_table --input ./summary_report/Summary.humann2.genefamilies.tsv --units cpm --output ./summary_report/Summary.humann2.genefamilies_cpm.tsv\n".format(**config_dict))
        # 
        f.write("python {script_dir}/megahit-quest-multi-result.py ./megahit_result \n".format(**config_dict))
        f.write("python {script_dir}/prokka_multi_result.py ./prokka_annotation \n".format(**config_dict))
        f.write("python {script_dir}/results2html.py \n".format(**config_dict))

    print("all finished!")

if __name__ == '__main__':
    main()
