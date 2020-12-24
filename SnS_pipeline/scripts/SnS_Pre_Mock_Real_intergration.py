import os
import sys
import argparse
import configparser

def usage():
    '''
    Usage:
        python3 SnS_Pre_Mock_Real_intergration.py -c SnS_Pre_Mock_Real_intergration.Comp-CR13-CD8T.config  -p Comp-CR13-CD8T

    Updated:
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
        'mock_samples' : cf.get('data', 'Mock'),
        'real_samples' : cf.get('data', 'Real'),
        }

    if len(config_dict['pre_samples']) == 0:
        pre_sample_list = []
    else:
        pre_sample_list = config_dict['pre_samples'].replace(" ","").rstrip(",").split(",")

    if len(config_dict['mock_samples']) == 0:
        mock_sample_list = []
    else:
        mock_sample_list = config_dict['mock_samples'].replace(" ","").rstrip(",").split(",")

    if len(config_dict['real_samples']) == 0:
        real_sample_list = []
    else:
        real_sample_list = config_dict['real_samples'].replace(" ","").rstrip(",").split(",")

    all_sample_list = pre_sample_list + mock_sample_list + real_sample_list
    print("your input sample names are: {}".format(all_sample_list))
    l1 = len(pre_sample_list)
    l2 = len(mock_sample_list)
    l3 = len(real_sample_list)
    print("pre_sample : {}\t mock_sample : {}\t Real_sample : {}".format(l1, l2, l3))
    sample_index_dict = {}
    sample_index_dict.setdefault('pre',  [x for x in range(1, l1+1)])
    sample_index_dict.setdefault('mock', [x for x in range(l1+1, l1+l2+1)])
    sample_index_dict.setdefault('real', [x for x in range(l1+l2+1, l1+l2+l3+1)])
    print("your sample_index dict: {}".format(sample_index_dict))
    # print("sample_index_dict should be: {'pre': [1, 2], 'mock': [3, 4], 'real': [5, 6]}")

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
    for sample in real_sample_list:
        tcrid = big_dict[sample].keys()
        real_sample_union_id += tcrid
    real_sample_union_id = set(real_sample_union_id)
    # print(real_sample_union_id)
    # print header first. 
    header = []
    header.append("TCR_id(Real_samples_union)")
    for each in pre_sample_list:
        header.append(each+"_count(Pre)")
        header.append(each+"_freq(Pre)")
    for each in mock_sample_list:
        header.append(each+"_count(Mock)")
        header.append(each+"_freq(Mock)")
    for each in real_sample_list:
        header.append(each+"_count(Real)")
        header.append(each+"_freq(Real)")
    for i in mock_sample_list:
        for j in pre_sample_list:
            header.append("Mock/Pre(" + i + "/" + j + ")")
    for i in real_sample_list:
        for j in pre_sample_list:
            header.append("Real/Pre(" + i + "/" + j + ")")
    for i in real_sample_list:
        for j in mock_sample_list:
            header.append("Real/Mock("+ i + "/" + j + ")")
    output.write("{}\n".format(",".join(header)))
    # print("finished write header")
    ##########################
    for each_id in real_sample_union_id:
        print_content = []
        print_content.append(each_id)
        for pre_sample in pre_sample_list:
            print_content.append(big_dict[pre_sample].get(each_id, {'perfect_count': '0'})['perfect_count'])
            print_content.append(big_dict[pre_sample].get(each_id, {'perfect_freq' : min_freq})['perfect_freq'])
        for mock_sample in mock_sample_list:
            print_content.append(big_dict[mock_sample].get(each_id, {'perfect_count': '0'})['perfect_count'])
            print_content.append(big_dict[mock_sample].get(each_id, {'perfect_freq' : min_freq})['perfect_freq'])
        for real_sample in real_sample_list:
            print_content.append(big_dict[real_sample].get(each_id, {'perfect_count': '0'})['perfect_count'])
            print_content.append(big_dict[real_sample].get(each_id, {'perfect_freq' : min_freq})['perfect_freq'])
        # print(print_content)

        # calculate mock/pre:
        foldchange_content = []
        for i in sample_index_dict['mock']:
            # print("mock index begin:{}".format(i))
            for j in sample_index_dict['pre']:
                # FC_value = calculate_division(print_content[i*2], print_content[j*2])
                # foldchange_content.append(FC_value)
                foldchange_content.append(str(float(print_content[i*2])/float(print_content[j*2])))
        for i in sample_index_dict['real']:
            for j in sample_index_dict['pre']:
                # FC_value = calculate_division(print_content[i*2], print_content[j*2])
                # foldchange_content.append(FC_value)
                foldchange_content.append(str(float(print_content[i*2])/float(print_content[j*2])))
        for i in sample_index_dict['real']:
            for j in sample_index_dict['mock']:
                # FC_value = calculate_division(print_content[i*2], print_content[j*2])
                # foldchange_content.append(FC_value)
                foldchange_content.append(str(float(print_content[i*2])/float(print_content[j*2])))
        output.write("{},{}\n".format(",".join(print_content), ",".join(foldchange_content)))
    output.close()
    print("all done!")

def calculate_division(a,b):
    a = float(a)
    b = float(b)
    if a == 0 and b == 0:
        v = str('0')
    elif a == 0 and b != 0:
        v = str(0.0000001/b)
    elif a != 0 and b == 0:
        v = str(a/0.0000001)
    elif a != 0 and b !=0:
        v = str(a/b)
    return v


if __name__ == '__main__':
    main()
