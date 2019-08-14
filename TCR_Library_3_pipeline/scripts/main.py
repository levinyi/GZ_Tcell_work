import sys
import os
import argparse
import ConfigParser


def _argparse():
    parser = argparse.ArgumentParser(description="This is description")
    parser.add_argument('-c', '--config', action='store', dest='config', default='config.txt', help="this is config file")
    parser.add_argument('-l', '--list', action='store', dest='data', default='data.list', help="this is data list file")
    return parser.parse_args()


def make_dir(*dir):
    for each in dir:
        if not os.path.exists(each):
            os.mkdir(each, 0755)


def addtwodimdict(thedict, key_a, key_b, val):
    ''' this is a function to add two dimetion dict '''
    if key_a in thedict:
        thedict[key_a].update({key_b: val})
    else:
        thedict.update({key_a: {key_b: val}})
    return thedict


def deal_rawdata(datalist, data_dir):
    '''
    read and deal with data.list file,
    '''
    rawdata_dict = {}
    with open(datalist) as f:
        for line in f:
            line = line.rstrip("\n")
            if line.startswith(("#", "\n")):
                continue
            sample_name, data = line.split("=")  # G13E3L12=/cygene/data/20190717_SH_G13_G23_G31/G13E3L12_S1_L001_R1_001.fastq.gz;/cygene/data/20190717_SH_G13_G23_G31/G13E3L12_S1_L001_R2_001.fastq.gz
            sample_name = sample_name.strip()  # if it contains space or something. then strip.
            R1, R2 = data.split(";")
            R1 = R1.strip()
            # check raw input read file exists in data.list.
            if not os.path.exists(R1):
                os.exit("Error: %s is not correct! please check your data.list " % R1)
            if not os.path.exists(R2):
                os.exit("Error: %s is not correct! please check your data.list " % R2)

            addtwodimdict(rawdata_dict, sample_name, R1, R2)

            # make a link to data dir.
            if not os.path.exists(data_dir + '/' + sample_name + '_R1.fq.gz'):
                os.system("ln -s %(R1)s %(data_dir)s/%(sample_name)s_R1.fq.gz" % {'R1': R1, 'data_dir': data_dir, 'sample_name': sample_name})
            if not os.path.exists(data_dir + '/' + sample_name + '_R2.fq.gz'):
                os.system("ln -s %(R2)s %(data_dir)s/%(sample_name)s_R2.fq.gz" % {'R2': R2, 'data_dir': data_dir, 'sample_name': sample_name})
    return rawdata_dict


def main():
    parser = _argparse()
    cf = ConfigParser.ConfigParser(allow_no_value=True)
    cf.read(parser.config)

    # check config file
    if not cf.has_section('Config'):
        os.exit("Error: your config file is not correct.")

    # read config file
    adapter1 = cf.get('Config', 'adapter1')
    adapter2 = cf.get('Config', 'adapter2')
    database = cf.get('Config', 'database')
    scripts = cf.get('Config', 'scripts')
    project_name = cf.get("Config", "project_name")
    merge_analysis = cf.get("Config", "merge_analysis")
    library_type = cf.get("Config", "library_type")

    if library_type == 'A_B':
        AB_ship = {'TRtype_R1': 'TRA', 'TRtype_R2': 'TRB', 'R_TRA': 'R1', 'R_TRB': 'R2'}
    elif library_type == 'B_A':
        AB_ship = {'TRtype_R1': 'TRB', 'TRtype_R2': 'TRA', 'R_TRA': 'R2', 'R_TRB': 'R1'}
    else:
        sys.exit("Error: your library_type is not correct, alternative option: 'A_B' or 'B_A'. ")

    # generate scripts
    # step 1 cutadapt
    work_path = os.path.abspath(".")

    # make all directories.
    project_dir = os.path.abspath(".") + '/' + project_name
    result_dir = project_dir + '/result'
    shell_dir = project_dir + '/shell'
    data_dir = project_dir + '/data'

    make_dir(project_dir, result_dir, shell_dir, data_dir)
    print("# Create work directory")

    # generate shell
    shell_name = shell_dir + '/work.' + project_name + '.sh'
    # only open a file so use try:finally to close.
    rawdata_dict = deal_rawdata(parser.data, data_dir)

    R1_dict = {}
    R2_dict = {}
    R1_list = []
    R2_list = []
    sample_list = []
    step1_shell = open(shell_name, "w")
    try:
        step1_shell.write("# cutadapt all samples\n")
        for sample_name, v in rawdata_dict.items():  # {G13E3L4 : {/path/to/sample_R1.fq.gz : /path/to/sample_R2.fq.gz}}
            sample_list.append(sample_name)
            for R1, R2 in v.items():
                R1_p1 = result_dir + '/' + sample_name + '_R1.p1.fq'
                R1_list.append(R1_p1)
                R1_dict[sample_name + '_R1'] = R1_p1
                R1_p2 = result_dir + '/' + sample_name + '_R1.p2.fq'

                R2_p1 = result_dir + '/' + sample_name + '_R2.p1.fq'
                R2_dict[sample_name + '_R2'] = R2
                R2_p2 = result_dir + '/' + sample_name + '_R2.p2.fq'
                R2_list.append(R2_p2)

                R1_logfile = result_dir + '/' + sample_name + '_R1.cutadapt.log'
                R2_logfile = result_dir + '/' + sample_name + '_R2.cutadapt.log'
                R1_statsfile = result_dir + '/' + sample_name + '_R1.cutadapt.stats'
                R2_statsfile = result_dir + '/' + sample_name + '_R2.cutadapt.stats'

                step1_shell.write("cutadapt -a adapter=%(adapter1)s -O 10 -o %(R1_p1)s -r %(R1_p2)s --info-file=%(R1_logfile)s %(R1_raw_data)s > %(R1_statsfile)s\n" % {"adapter1": adapter1, 'R1_p1': R1_p1, 'R1_p2': R1_p2, 'R1_logfile': R1_logfile, 'R1_raw_data': R1, 'R1_statsfile': R1_statsfile})
                step1_shell.write("cutadapt -g adapter=%(adapter2)s -O 10 -o %(R2_p1)s -r %(R2_p2)s --info-file=%(R2_logfile)s %(R2_raw_data)s > %(R2_statsfile)s\n" % {"adapter2": adapter2, 'R2_p1': R2_p1, 'R2_p2': R2_p2, 'R2_logfile': R2_logfile, 'R2_raw_data': R2, 'R2_statsfile': R2_statsfile})

                step1_shell.write("python  %(scripts)s/deal_cutadapt_log.py -l %(R1_logfile)s -d %(result_dir)s  \n" % {'scripts': scripts, 'result_dir': result_dir, 'R1_logfile': R1_logfile})
                step1_shell.write("python  %(scripts)s/deal_cutadapt_log.py -l %(R2_logfile)s -d %(result_dir)s  \n" % {'scripts': scripts, 'result_dir': result_dir, 'R2_logfile': R2_logfile})

        # check merge
        if merge_analysis == 'True':
            print "# merge analysis : True"
            all_R1 = " ".join(R1_list)
            all_R2 = " ".join(R2_list)
            step1_shell.write("\n# merge analysis : True\n")
            step1_shell.write("cat %(all_R1)s >%(result_dir)s/%(project_name)s_R1.p1.fq\n" % {'all_R1': all_R1, 'result_dir': result_dir, 'project_name': project_name})
            step1_shell.write("cat %(all_R2)s >%(result_dir)s/%(project_name)s_R2.p2.fq\n" % {'all_R2': all_R2, 'result_dir': result_dir, 'project_name': project_name})

            # store a dict for next step.
            new_data_dict = {}
            new_data_dict[project_name + '_R1'] = "%(result_dir)s/%(project_name)s_R1.p1.fq" % {'result_dir': result_dir, 'project_name': project_name}
            new_data_dict[project_name + '_R2'] = "%(result_dir)s/%(project_name)s_R2.p2.fq" % {'result_dir': result_dir, 'project_name': project_name}
        else:
            print "# merg analysis : False"
            step1_shell.write("\n# merge analysis : False\n")
            new_data_dict = dict(R1_dict, **R2_dict)

        # step 2  mixcr
        step1_shell.write("\n# mixcr\n")
        for prefix, new_data in new_data_dict.items():
            make_dir(result_dir + '/' + prefix)
            step1_shell.write("mixcr analyze shotgun --align \"-OsaveOriginalReads=true\" --species hs --starting-material rna --receptor-type tcr -r %(result_dir)s/%(prefix)s.report %(fastq)s %(result_dir)s/%(prefix)s.mixcr.out \n" % {'prefix': prefix, 'fastq': new_data, 'result_dir': result_dir})
            step1_shell.write("mixcr exportReadsForClones -s %(result_dir)s/%(prefix)s.mixcr.out.clna %(result_dir)s/%(prefix)s/%(prefix)s.\n\n" % {'result_dir': result_dir, 'prefix': prefix})

        # generate list file
        with open(result_dir + "/list", "w") as f:
            for each in sample_list:
                f.write("%s\n" % each)

        print("# Your library_type is %s" % library_type)
        step1_shell.write("# library_type : %s\n\n" % library_type)

        variable_dict = {'project_dir': project_dir, 'scripts': scripts, 'project_name': project_name, 'result_dir': result_dir, 'TRtype_R1': AB_ship['TRtype_R1'],
                         'TRtype_R2': AB_ship['TRtype_R2'], 'library_type': library_type}
        if merge_analysis == 'True':
            step1_shell.write("perl %(scripts)s/01_separate_reads.pl %(result_dir)s/list R1 %(TRtype_R1)s %(project_name)s \n" % variable_dict)
            step1_shell.write("perl %(scripts)s/01_separate_reads.pl %(result_dir)s/list R2 %(TRtype_R2)s %(project_name)s \n" % variable_dict)
            step1_shell.write("perl %(scripts)s/02_identify_pair.pl  %(result_dir)s/list %(library_type)s\n\n" % variable_dict)
            step1_shell.write("python %(scripts)s/statistic_basic_info.py %(project_dir)s > %(result_dir)s/%(project_name)s.basic_infomation.xls\n\n" % variable_dict)

            plot_dir = result_dir + '/plot'
            make_dir(plot_dir)

            step1_shell.write("# draw plots: \n\n")
            for sample in sample_list:
                freq = result_dir + '/' + sample + '.pairs.freq'
                r1_reads_file = result_dir + '/' + sample + '_R1.' + AB_ship['TRtype_R1'] + '.reads'
                r2_reads_file = result_dir + '/' + sample + '_R2.' + AB_ship['TRtype_R2'] + '.reads'

                variable_dict = {'freq': freq, 'scripts': scripts, 'result_dir': result_dir, 'sample': sample, 'project_name': project_name, 'plot_dir': plot_dir,
                                 'R_TRA': AB_ship['R_TRA'], 'R_TRB': AB_ship['R_TRB'], 'TRtype_R1': AB_ship['TRtype_R1'], 'TRtype_R2': AB_ship['TRtype_R2'], 'r1_reads_file': r1_reads_file, 'r2_reads_file': r2_reads_file}
                # reads 2 clonotype
                # print r1_reads_file,r2_reads_file
                step1_shell.write("python %(scripts)s/reads2relaxed_clonotype.py %(result_dir)s/%(project_name)s_R1.mixcr.out.clonotypes.%(TRtype_R1)s.txt %(r1_reads_file)s > %(result_dir)s/%(sample)s_R1.%(TRtype_R1)s.relaxed.clonotype.reads.xls\n" % variable_dict)
                step1_shell.write("python %(scripts)s/reads2relaxed_clonotype.py %(result_dir)s/%(project_name)s_R2.mixcr.out.clonotypes.%(TRtype_R2)s.txt %(r2_reads_file)s > %(result_dir)s/%(sample)s_R2.%(TRtype_R2)s.relaxed.clonotype.reads.xls\n" % variable_dict)

                # 01.matrix plot
                step1_shell.write("Rscript %(scripts)s/matrix_raster.pairing.R  %(plot_dir)s  %(freq)s \n" % variable_dict)
                # 02.top pair dominant # input:freq and dominant.txt
                step1_shell.write("python  %(scripts)s/top_pair_dominant.py %(freq)s >%(result_dir)s/%(sample)s.top_pair_dominant.txt \n" % variable_dict)
                # for 03 plot # input :sharpen.data
                step1_shell.write("perl  %(scripts)s/03_sharpen_pair.pl %(sample)s 3 0\n" % variable_dict)
                # generate M1 table
                step1_shell.write("python  %(scripts)s/identify_dominant_pairs_M1_for_new_method.py %(freq)s %(freq)s.M1.xls %(result_dir)s/%(project_name)s_%(R_TRA)s.mixcr.out.clonotypes.TRA.txt %(result_dir)s/%(project_name)s_%(R_TRB)s.mixcr.out.clonotypes.TRB.txt\n" % variable_dict)
                step1_shell.write("python  %(scripts)s/M1_relaxed.py %(freq)s.M1.xls > %(freq)s.relaxed.M1.xls\n" % variable_dict)
                # reads number of each gene
                step1_shell.write("python  %(scripts)s/gene_reads_number.py %(freq)s.M1.xls 1 >%(result_dir)s/%(sample)s.filter.M1.TRV.reads.number.xls\n" % variable_dict)
                step1_shell.write("python  %(scripts)s/gene_reads_number.py %(freq)s.M1.xls 2 >%(result_dir)s/%(sample)s.total.M1.TRV.reads.number.xls\n" % variable_dict)
                # reads length of after cutadapt
                step1_shell.write("python  %(scripts)s/statistic_reads_length.py %(result_dir)s/%(sample)s_R1.p1.fq >%(result_dir)s/%(sample)s_R1.p1.fq.length.txt\n" % variable_dict)
                step1_shell.write("python  %(scripts)s/statistic_reads_length.py %(result_dir)s/%(sample)s_R2.p2.fq >%(result_dir)s/%(sample)s_R2.p2.fq.length.txt\n\n" % variable_dict)

            # draw point plot of dominant pairs.
            step1_shell.write("Rscript %(scripts)s/dominant.point.R %(plot_dir)s %(result_dir)s/*.top_pair_dominant.txt\n" % variable_dict)
            # draw read fraction of dominant pairs.
            step1_shell.write("Rscript %(scripts)s/03.read_frac_for_dom_pairs.R %(plot_dir)s  %(result_dir)s/*.sharpen.data\n" % variable_dict)
            # draw bar plot of reads number.
            step1_shell.write("Rscript %(scripts)s/gene_reads_number_barplot.R  %(plot_dir)s  %(result_dir)s/*.M1.TRV.reads.number.xls\n" % variable_dict)
            # draw histogram of reads length.
            step1_shell.write("Rscript %(scripts)s/statistic_reads_length.R %(plot_dir)s %(result_dir)s/*.length.txt\n\n" % variable_dict)
            # print("Congratulations! All projects directories were created!")
            # print("Go %s and run the shell in shell directory " % project_name)
            step1_shell.close()
        else:  # merg analysis : False
            for prefix in new_data_dict:
                step1_shell.write("python %(scripts)s/slim.mixcr.py %(result_dir)s/%(prefix)s.mixcr.out.clonotypes.TRA.txt > %(result_dir)s/%(prefix)s.clonotypes.reads.xls \n" % {'scripts': scripts, 'result_dir': result_dir, 'prefix': prefix})
                step1_shell.write("python %(scripts)s/statistic_basic_info.py %(project_dir)s > %(result_dir)s/%(project_name)s.basic_infomation.xls\n\n" % variable_dict)
    finally:
        step1_shell.close()

    print("Congratulations! All projects directories were created!")
    print("Go %s and run the shell in shell directory " % project_name)

if __name__ == '__main__':
    main()
