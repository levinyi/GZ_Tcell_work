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
        'tra_fastqs': cf.get('Config', 'tra_fastqs'),
        'tra_fastq_rd1' : cf.get('Config', 'tra_fastqs').split()[0],
        'trb_fastqs': cf.get('Config', 'trb_fastqs'),
        'trb_fastq_rd1': cf.get('Config', 'trb_fastqs').split()[0],
        'tra_barcode_file': cf.get('Config', 'tra_barcode_file'),
        'trb_barcode_file': cf.get('Config', 'trb_barcode_file'),
        
        # project name
        'project_name' : cf.get("Config", "project_name"),
        'tra_sample_name'  : cf.get("Config", "tra_sample_name"),
        'trb_sample_name'  : cf.get("Config", "trb_sample_name"),
        'merged_sample_name' : cf.get("Config", "merged_sample_name"),
        'well_number': cf.get("Config", "well_number"),
        'expected_cell_number': cf.get("Config", "expected_cell_number"),

        # alignment
        'mixcr': cf.get('Config', "mixcr"),
        'scripts_dir': cf.get('Config', "scripts_dir"),
        # 'database_dir': cf.get('Config', "database_dir"),
    }
    #make directories:
    project_dir = os.path.abspath(".") + '/' + config_dict['project_name']
    make_dir(project_dir)
    print("# Create work directory")

    # generate shell
    shell_name = project_dir + '/work.' + config_dict['project_name'] + '.sh'
    # only open a file so use try:finally to close.

    with open(shell_name,"w") as f:
        # align
        f.write("{mixcr} -Xmx64G -Xms64G analyze shotgun --align -OsaveOriginalReads=true --species hs --starting-material rna --receptor-type tcr  -r {tra_sample_name}.report  {tra_fastqs} {tra_sample_name}.mixcr.out \n".format(**config_dict))
        f.write("{mixcr} -Xmx64G -Xms64G analyze shotgun --align -OsaveOriginalReads=true --species hs --starting-material rna --receptor-type tcr  -r {trb_sample_name}.report  {trb_fastqs} {trb_sample_name}.mixcr.out \n".format(**config_dict))
        f.write("mkdir -p {tra_sample_name} {trb_sample_name} \n".format(**config_dict))
        f.write("{mixcr} exportReadsForClones -s {tra_sample_name}.mixcr.out.clna {tra_sample_name}/{tra_sample_name}.\n".format(**config_dict))
        f.write("{mixcr} exportReadsForClones -s {trb_sample_name}.mixcr.out.clna {trb_sample_name}/{trb_sample_name}.\n".format(**config_dict))

        f.write("python3 {scripts_dir}/0.read2barcode.py {tra_barcode_file} {tra_fastq_rd1} {tra_sample_name}\n".format(**config_dict))
        f.write("python3 {scripts_dir}/0.read2barcode.py {trb_barcode_file} {trb_fastq_rd1} {trb_sample_name}\n".format(**config_dict))
        f.write("less {tra_sample_name}.barcode.total.rate.xls |head -{well_number} | awk '{{print$1\"\\t\"$2}}' >{tra_sample_name}.barcode.txt \n".format(**config_dict))
        f.write("less {trb_sample_name}.barcode.total.rate.xls |head -{well_number} | awk '{{print$1\"\\t\"$2}}' >{trb_sample_name}.barcode.txt \n".format(**config_dict))
        f.write("mv {tra_sample_name}.barcode.total.rate.xls {tra_sample_name}.barcode.total.rate.raw.xls \n".format(**config_dict))
        f.write("mv {trb_sample_name}.barcode.total.rate.xls {trb_sample_name}.barcode.total.rate.raw.xls \n".format(**config_dict))
        f.write("python3 {scripts_dir}/0.read2barcode.py {tra_sample_name}.barcode.txt {tra_fastq_rd1} {tra_sample_name}\n".format(**config_dict))
        f.write("python3 {scripts_dir}/0.read2barcode.py {trb_sample_name}.barcode.txt {trb_fastq_rd1} {trb_sample_name}\n".format(**config_dict))

        f.write("python {scripts_dir}/1.mixcr2barcode.py {tra_sample_name}.barcode.read.txt {tra_sample_name}.mixcr.out.clonotypes.TRA.txt TRA\n".format(**config_dict))
        f.write("python {scripts_dir}/1.mixcr2barcode.py {trb_sample_name}.barcode.read.txt {trb_sample_name}.mixcr.out.clonotypes.TRB.txt TRB\n".format(**config_dict))
        f.write("cat {tra_sample_name}.TRA.clone.reads.barcode.txt {trb_sample_name}.TRB.clone.reads.barcode.txt  >Total.{merged_sample_name}.TRAB.clone.reads.barcode.txt \n".format(**config_dict))

        ##########################
        f.write("python3 {scripts_dir}/mixcr_filter_clonotype_counts.py {tra_sample_name}.mixcr.out.clonotypes.TRA.txt > {tra_sample_name}.mixcr.out.clonotypes.TRA.count.txt \n".format(**config_dict))
        f.write("python3 {scripts_dir}/mixcr_filter_clonotype_counts.py {trb_sample_name}.mixcr.out.clonotypes.TRB.txt > {trb_sample_name}.mixcr.out.clonotypes.TRB.count.txt \n".format(**config_dict))
        
        f.write("python3 {scripts_dir}/2.transTalbe2matrix.py -i Total.{merged_sample_name}.TRAB.clone.reads.barcode.txt --cloneCountThreshold 19 --overallThreshold 5 --column_median 10 --row_median 10 --filter_noise_wells \n".format(**config_dict))
        # f.write("python3 {scripts_dir}/2.1.optional_filter_median.py Total.{merged_sample_name}.TRAB.clone.reads.barcode.txt \n".format(**config_dict))
        f.write("Rscript {scripts_dir}/2.2.draw.readCount.plot.R *.csv \n".format(**config_dict))

        f.write("python3 {scripts_dir}/5.permutation_test.multi_process.py Total.{merged_sample_name}.TRAB.wells.boole.matrix.Filtered.csv  {merged_sample_name}.permutation_test.10000.out.Filtered.txt \n".format(**config_dict))
        f.write("mkdir -p {merged_sample_name}_Shotgun/FromA2B {merged_sample_name}_Shotgun/FromB2A \n\n".format(**config_dict))
        ######################### FromA2B
        f.write("python3 {scripts_dir}/3.calculate.p-value.FromA2B.py Total.{merged_sample_name}.TRAB.clone.reads.barcode.Filtered.txt 50 >{merged_sample_name}_Shotgun/FromA2B/Total.pairs.FromA2B.threshold.50.xls \n".format(**config_dict))
        f.write("python3 {scripts_dir}/3.calculate.p-value.FromA2B.py Total.{merged_sample_name}.TRAB.clone.reads.barcode.Filtered.txt 10 >{merged_sample_name}_Shotgun/FromA2B/Total.pairs.FromA2B.threshold.10.xls \n".format(**config_dict))
        f.write("python3 {scripts_dir}/6.permutation_true.v3.py {merged_sample_name}.permutation_test.10000.out.Filtered.txt {merged_sample_name}_Shotgun/FromA2B/Total.pairs.FromA2B.threshold.50.xls {expected_cell_number} > {merged_sample_name}_Shotgun/FromA2B/Total.pairs.FromA2B.threshold.50.add.xls \n".format(**config_dict))
        f.write("python3 {scripts_dir}/6.permutation_true.v3.py {merged_sample_name}.permutation_test.10000.out.Filtered.txt {merged_sample_name}_Shotgun/FromA2B/Total.pairs.FromA2B.threshold.10.xls {expected_cell_number} > {merged_sample_name}_Shotgun/FromA2B/Total.pairs.FromA2B.threshold.10.add.xls \n".format(**config_dict))
        f.write("less {merged_sample_name}_Shotgun/FromA2B/Total.pairs.FromA2B.threshold.50.add.xls |awk '$7>=5 && $10<0.1' >{merged_sample_name}_Shotgun/FromA2B/Total.pairs.FromA2B.threshold.50.add.filtered.xls \n".format(**config_dict))
        f.write("less {merged_sample_name}_Shotgun/FromA2B/Total.pairs.FromA2B.threshold.10.add.xls |awk '$7>=5 && $10<0.1' >{merged_sample_name}_Shotgun/FromA2B/Total.pairs.FromA2B.threshold.10.add.filtered.xls \n".format(**config_dict))
        f.write("less {merged_sample_name}_Shotgun/FromA2B/Total.pairs.FromA2B.threshold.50.add.xls |awk '$7>=3 && $10<0.1' >{merged_sample_name}_Shotgun/FromA2B/Total.pairs.FromA2B.threshold.50.add.filtered.sharedwells.3.xls \n".format(**config_dict))
        f.write("less {merged_sample_name}_Shotgun/FromA2B/Total.pairs.FromA2B.threshold.10.add.xls |awk '$7>=3 && $10<0.1' >{merged_sample_name}_Shotgun/FromA2B/Total.pairs.FromA2B.threshold.10.add.filtered.sharedwells.3.xls \n\n".format(**config_dict))
        
        ########################## FromB2A
        f.write("python3 {scripts_dir}/3.calculate.p-value.FromB2A.py Total.{merged_sample_name}.TRAB.clone.reads.barcode.Filtered.txt 50 >{merged_sample_name}_Shotgun/FromB2A/Total.pairs.FromB2A.threshold.50.xls \n".format(**config_dict))
        f.write("python3 {scripts_dir}/3.calculate.p-value.FromB2A.py Total.{merged_sample_name}.TRAB.clone.reads.barcode.Filtered.txt 10 >{merged_sample_name}_Shotgun/FromB2A/Total.pairs.FromB2A.threshold.10.xls \n".format(**config_dict))
        f.write("python3 {scripts_dir}/6.permutation_true.v3.py {merged_sample_name}.permutation_test.10000.out.Filtered.txt {merged_sample_name}_Shotgun/FromB2A/Total.pairs.FromB2A.threshold.50.xls {expected_cell_number} > {merged_sample_name}_Shotgun/FromB2A/Total.pairs.FromB2A.threshold.50.add.xls \n".format(**config_dict))
        f.write("python3 {scripts_dir}/6.permutation_true.v3.py {merged_sample_name}.permutation_test.10000.out.Filtered.txt {merged_sample_name}_Shotgun/FromB2A/Total.pairs.FromB2A.threshold.10.xls {expected_cell_number} > {merged_sample_name}_Shotgun/FromB2A/Total.pairs.FromB2A.threshold.10.add.xls \n".format(**config_dict))
        f.write("less {merged_sample_name}_Shotgun/FromB2A/Total.pairs.FromB2A.threshold.50.add.xls |awk '$7>=5 && $10<0.1' >{merged_sample_name}_Shotgun/FromB2A/Total.pairs.FromB2A.threshold.50.add.filtered.xls \n".format(**config_dict))
        f.write("less {merged_sample_name}_Shotgun/FromB2A/Total.pairs.FromB2A.threshold.10.add.xls |awk '$7>=5 && $10<0.1' >{merged_sample_name}_Shotgun/FromB2A/Total.pairs.FromB2A.threshold.10.add.filtered.xls \n".format(**config_dict))
        f.write("less {merged_sample_name}_Shotgun/FromB2A/Total.pairs.FromB2A.threshold.50.add.xls |awk '$7>=3 && $10<0.1' >{merged_sample_name}_Shotgun/FromB2A/Total.pairs.FromB2A.threshold.50.add.filtered.sharedwells.3.xls \n".format(**config_dict))
        f.write("less {merged_sample_name}_Shotgun/FromB2A/Total.pairs.FromB2A.threshold.10.add.xls |awk '$7>=3 && $10<0.1' >{merged_sample_name}_Shotgun/FromB2A/Total.pairs.FromB2A.threshold.10.add.filtered.sharedwells.3.xls \n\n".format(**config_dict))
        
        ########################### sub
        f.write("python3 {scripts_dir}/9.compare.A2B.v3.py {merged_sample_name}_Shotgun/FromA2B/Total.pairs.FromA2B.threshold.10.add.filtered.sharedwells.3.xls {merged_sample_name}_Shotgun/FromB2A/Total.pairs.FromB2A.threshold.10.add.filtered.sharedwells.3.xls > {merged_sample_name}_Shotgun/A2B_B2A_detail.threshold.10.sharedwells.3.xls \n".format(**config_dict))
        f.write("python3 {scripts_dir}/9.compare.A2B.v3.py {merged_sample_name}_Shotgun/FromA2B/Total.pairs.FromA2B.threshold.50.add.filtered.sharedwells.3.xls {merged_sample_name}_Shotgun/FromB2A/Total.pairs.FromB2A.threshold.50.add.filtered.sharedwells.3.xls > {merged_sample_name}_Shotgun/A2B_B2A_detail.threshold.50.sharedwells.3.xls \n".format(**config_dict))
        f.write("python3 {scripts_dir}/10.merge.A2B.B2A.v1.py {merged_sample_name}_Shotgun/FromA2B/Total.pairs.FromA2B.threshold.10.add.filtered.sharedwells.3.xls {merged_sample_name}_Shotgun/FromB2A/Total.pairs.FromB2A.threshold.10.add.filtered.sharedwells.3.xls  {merged_sample_name}_Shotgun/A2B_B2A_detail.threshold.10.sharedwells.3.new.xlsx \n".format(**config_dict))
        f.write("python3 {scripts_dir}/10.merge.A2B.B2A.v1.py {merged_sample_name}_Shotgun/FromA2B/Total.pairs.FromA2B.threshold.50.add.filtered.sharedwells.3.xls {merged_sample_name}_Shotgun/FromB2A/Total.pairs.FromB2A.threshold.50.add.filtered.sharedwells.3.xls  {merged_sample_name}_Shotgun/A2B_B2A_detail.threshold.50.sharedwells.3.new.xlsx \n".format(**config_dict))

    print("all finished!")


if __name__ == '__main__':
    main()
