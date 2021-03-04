import os
import sys
import argparse
import configparser
import numpy as np

def usage():
    '''
    Usage:
        python3 SnS_Pre_Mock_Real_intergration.py -c SnS_Pre_Mock_Real_intergration.Comp-CR13-CD8T.config  -p Comp-CR13-CD8T

    Updated:
    20210301: updated a new version, allow multiple mock sample(test for 6 mock sample),and multiple real sample.
                Using specific samples of average frequency for next fold change calculation.
    20201222: fix a bug. when config file add space between two samples.
    '''

def _argparse():
    parser = argparse.ArgumentParser(description="This is description")
    parser.add_argument('-c', '--config',  action='store', dest='config', default='config.txt', help="this is config file")
    parser.add_argument('-p', '--project', action='store', dest='project', default='my_project', help='this is project')
    return parser.parse_args()


def addtwodimdict(thedict, key_a, key_b, val):
    ''' this is a function to add two dimetion dict '''
    if key_a in thedict:
        thedict[key_a].update({key_b: val})
    else:
        thedict.update({key_a: {key_b: val}})
    return thedict


def read_freq_file(freq_file, min_freq):
    adict = {}
    with open(freq_file, "r") as f:
        for line in f:
            if line.startswith("pair"):
                continue
            line = line.rstrip("\n")
            pair, cdr3acc,cdr3freq, perfect_count, perfect_freq = line.split()
            if float(perfect_freq) == 0:
                addtwodimdict(adict, pair, 'perfect_freq', min_freq)
            else:
                addtwodimdict(adict, pair, 'perfect_freq', perfect_freq)
            addtwodimdict(adict, pair, 'perfect_count', perfect_count)
    return adict

def calculate_average_freq(sample_index_dict, x, print_content, min_freq):
    freq_list = []
    for i in sample_index_dict[x]:
        freq = print_content[i*2]
        if freq != min_freq:
            freq_list.append(float(freq))
    if len(freq_list) == 0:
        average_freq = min_freq
    else:
        average_freq = str(float(np.mean(freq_list)))
    return average_freq


def main():
    parser = _argparse()
    cf = configparser.ConfigParser(allow_no_value=False, )
    cf.read(parser.config)
    if not cf.has_section('data'):
        os.exit("Error: your config file is not correct.")
    project_name = parser.project
    output = open(project_name + '.csv', "w")

    config_dict = {
        'pre_samples' : cf.get('data','Pre'),
        'mock1_samples' : cf.get('data', 'Mock-1'),
        'mock2_samples' : cf.get('data', 'Mock-2'),
        'real1_samples' : cf.get('data', 'Real-1'),
        'real2_samples' : cf.get('data', 'Real-2'),
        }

    if len(config_dict['pre_samples']) == 0:
        pre_sample_list = []
    else:
        pre_sample_list = config_dict['pre_samples'].replace(" ","").rstrip(",").split(",")
    if len(config_dict['mock1_samples']) == 0:
        mock1_sample_list = []
    else:
        mock1_sample_list = config_dict['mock1_samples'].replace(" ","").rstrip(",").split(",")
    if len(config_dict['mock2_samples']) == 0:
        mock2_sample_list = []
    else:
        mock2_sample_list = config_dict['mock2_samples'].replace(" ","").rstrip(",").split(",")
    if len(config_dict['real1_samples']) == 0:
        real1_sample_list = []
    else:
        real1_sample_list = config_dict['real1_samples'].replace(" ","").rstrip(",").split(",")
    if len(config_dict['real2_samples']) == 0:
        real2_sample_list = []
    else:
        real2_sample_list = config_dict['real2_samples'].replace(" ","").rstrip(",").split(",")

    all_sample_list = pre_sample_list + mock1_sample_list + mock2_sample_list + real1_sample_list + real2_sample_list
    print("your input sample names are: {}".format(all_sample_list))
###########################################################
    min_freq = '0.0000001'
    files = os.listdir(path='.')
    big_dict = {}
    for sample in all_sample_list:
        sample_file = [each for each in files if each.startswith(sample+'.pair.acc_freq')]
        print("detected file: {}".format(sample_file))
        if len(sample_file) > 1:
            sys.exit("Error!")
        try:
            sample_file = sample_file[0]
        except IndexError:
            sys.exit("Error: Wrong sample name in your config file: {}. Could not find {}.pair.acc_freq*.txt".format(sample,sample))
        adict = read_freq_file(sample_file, min_freq)
        big_dict[sample] = adict
    # print(big_dict)
    real_sample_union_id = []
    for sample in real1_sample_list+real2_sample_list:
        tcrid = big_dict[sample].keys()
        real_sample_union_id += tcrid
    real_sample_union_id = set(real_sample_union_id)
    # print(real_sample_union_id)
    ####################################################################
    # print header first. 
    header = []
    header.append("TCR_id(Real_samples_union)")
    for each in pre_sample_list:
        header.append(each+"_count(Pre)")
        header.append(each+"_freq(Pre)")
    for each in mock1_sample_list:
        header.append(each+"_count(Mock1)")
        header.append(each+"_freq(Mock1)")
    for each in mock2_sample_list:
        header.append(each+"_count(Mock2)")
        header.append(each+"_freq(Mock2)")
    for each in real1_sample_list:
        header.append(each+"_count(Real1)")
        header.append(each+"_freq(Real1)")
    for each in real2_sample_list:
        header.append(each+"_count(Real2)")
        header.append(each+"_freq(Real2)")
    
    header.append("mock1_freq_average({})".format(";".join(mock1_sample_list)))
    header.append("mock2_freq_average({})".format(";".join(mock2_sample_list)))
    header.append("real1_freq_average({})".format(";".join(real1_sample_list)))
    header.append("real2_freq_average({})".format(";".join(real2_sample_list)))
    new_mock_sample_list = ['mock1_ave', 'mock2_ave']
    new_real_sample_list = ['real1_ave', 'real2_ave']
    for i in new_mock_sample_list:
        for j in pre_sample_list:
            header.append("Mock/Pre(" + i + "/" + j + ")")
    for i in new_real_sample_list:
        for j in pre_sample_list:
            header.append("Real/Pre(" + i + "/" + j + ")")
    for i in new_real_sample_list:
        for j in new_mock_sample_list:
            header.append("Real/Mock("+ i + "/" + j + ")")
    output.write("{}\n".format(",".join(header)))
    # print("finished write header")

    ####################################################
    l1 = len(pre_sample_list)
    l2 = len(mock1_sample_list)
    l3 = len(mock2_sample_list)
    l4 = len(real1_sample_list)
    l5 = len(real2_sample_list)
    print("pre_sample : {}\t mock1_sample : {}\t mock2_sample : {}\t Real1_sample : {} \tReal2_sample : {}".format(l1, l2, l3,l4,l5))
    sample_index_dict = {}
    sample_index_dict.setdefault('pre',  [x for x in range(1, l1+1)])
    sample_index_dict.setdefault('mock1', [x for x in range(l1+1, l1+l2+1)])
    sample_index_dict.setdefault('mock2', [x for x in range(l1+l2+1, l1+l2+l3+1)])
    sample_index_dict.setdefault('real1', [x for x in range(l1+l2+l3+1, l1+l2+l3+l4+1)])
    sample_index_dict.setdefault('real2', [x for x in range(l1+l2+l3+l4+1, l1+l2+l3+l4+l5+1)])
    print("your sample_index dict: {}".format(sample_index_dict))
    print("sample_index_dict should be: {'pre': [1, 2], 'mock1':[3,4,5], 'mock2':[6,7,8], 'real1':[9,10,11], 'real2':[12,13,14]}")
    
    ####################################################
    for each_id in real_sample_union_id:
        print_content = []
        print_content.append(each_id)
        for pre_sample in pre_sample_list:
            print_content.append(big_dict[pre_sample].get(each_id, {'perfect_count': '0'})['perfect_count'])
            print_content.append(big_dict[pre_sample].get(each_id, {'perfect_freq' : min_freq})['perfect_freq'])
        for mock_sample in mock1_sample_list:
            print_content.append(big_dict[mock_sample].get(each_id, {'perfect_count': '0'})['perfect_count'])
            print_content.append(big_dict[mock_sample].get(each_id, {'perfect_freq' : min_freq})['perfect_freq'])
        for mock_sample in mock2_sample_list:
            print_content.append(big_dict[mock_sample].get(each_id, {'perfect_count': '0'})['perfect_count'])
            print_content.append(big_dict[mock_sample].get(each_id, {'perfect_freq' : min_freq})['perfect_freq'])
        for real_sample in real1_sample_list:
            print_content.append(big_dict[real_sample].get(each_id, {'perfect_count': '0'})['perfect_count'])
            print_content.append(big_dict[real_sample].get(each_id, {'perfect_freq' : min_freq})['perfect_freq'])
        for real_sample in real2_sample_list:
            print_content.append(big_dict[real_sample].get(each_id, {'perfect_count': '0'})['perfect_count'])
            print_content.append(big_dict[real_sample].get(each_id, {'perfect_freq' : min_freq})['perfect_freq'])
        # print(print_content) 

        # calculate mock1 average:
        mock1_ave_freq = calculate_average_freq(sample_index_dict,'mock1', print_content, min_freq)
        mock2_ave_freq = calculate_average_freq(sample_index_dict,'mock2', print_content, min_freq)
        real1_ave_freq = calculate_average_freq(sample_index_dict,'real1', print_content, min_freq)
        real2_ave_freq = calculate_average_freq(sample_index_dict,'real2', print_content, min_freq)
        print_content.append(mock1_ave_freq)
        print_content.append(mock2_ave_freq)
        print_content.append(real1_ave_freq)
        print_content.append(real2_ave_freq)
        print_content_len = len(print_content)
        # print(print_content_len)
        new_sample_index_dict = {}
        new_sample_index_dict.setdefault('pre',  [x for x in range(1, l1+1)])
        new_sample_index_dict.setdefault('mock', [x for x in range(print_content_len-4, print_content_len-4+2)])
        new_sample_index_dict.setdefault('real', [x for x in range(print_content_len-4+2, print_content_len-4+4)])
        # print("new_sample_index: {}".format(new_sample_index_dict))
        # print("sample_index_dict should be: {'pre': [1, 2], 'mock':[29,30], 'real':[31,32]}")
        # calculate mock/pre:
        foldchange_content = []
        for i in new_sample_index_dict['mock']:
            for j in new_sample_index_dict['pre']:
                foldchange_content.append(str(float(print_content[i])/float(print_content[j*2])))
        for i in new_sample_index_dict['real']:
            for j in new_sample_index_dict['pre']:
                foldchange_content.append(str(float(print_content[i])/float(print_content[j*2])))
        for i in new_sample_index_dict['real']:
            for j in new_sample_index_dict['mock']:
                foldchange_content.append(str(float(print_content[i])/float(print_content[j])))
        output.write("{},{}\n".format(",".join(print_content), ",".join(foldchange_content)))
    output.close()
    print("all done!")


if __name__ == '__main__':
    main()
