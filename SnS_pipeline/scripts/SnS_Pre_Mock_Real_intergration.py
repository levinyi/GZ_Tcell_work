import os
import sys
import argparse
import configparser


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


def read_freq_file(freq_file):
    adict = {}
    with open(freq_file, "r") as f:
        for line in f:
            if line.startswith("pair"):
                continue
            line = line.rstrip("\n")
            pair, cdr3acc,cdr3freq, perfect_count, perfect_freq = line.split()
            addtwodimdict(adict, pair, 'perfect_count', perfect_count)
            addtwodimdict(adict, pair, 'perfect_freq', perfect_freq)
    return adict


def main():
    parser = _argparse()
    cf = configparser.ConfigParser(allow_no_value=True)
    cf.read(parser.config)
    if not cf.has_section('data'):
        os.exit("Error: your config file is not correct.")
    project_name = parser.project

    config_dict = {
        'pre_samples' : cf.get('data','Pre'),
        'mock_samples' : cf.get('data', 'Mock'),
        'real_samples' : cf.get('data', 'Real'),
        }
    
    pre_sample_list = config_dict['pre_samples'].split(",")
    mock_sample_list = config_dict['mock_samples'].split(",")
    real_sample_list = config_dict['real_samples'].split(",")
    all_sample_list = pre_sample_list + mock_sample_list + real_sample_list
    l1 = len(pre_sample_list)
    l2 = len(mock_sample_list)
    l3 = len(real_sample_list)
    
    sample_index_dict = {}
    sample_index_dict.setdefault('pre',  [x for x in range(1, l1+1)])
    sample_index_dict.setdefault('mock', [x for x in range(l1+1, l1+l2+1)])
    sample_index_dict.setdefault('real', [x for x in range(l1+l2+1, l1+l2+l3+1)])
    # print(sample_index_dict)
    # sample_index_dict: {'pre': [1, 2], 'mock': [3, 4], 'real': [5, 6]}

    files = os.listdir(path='.')
    big_dict = {}
    for sample in all_sample_list:
        sample_file = [each for each in files if each.startswith(sample+'.pair.acc_freq')]
        if len(sample_file) > 1:
            sys.exit("Error!")
        sample_file = sample_file[0]
        adict = read_freq_file(sample_file)
        big_dict[sample] = adict
    
    real_sample_union_id = []
    for sample in real_sample_list:
        tcrid = big_dict[sample].keys()
        real_sample_union_id += tcrid
    real_sample_union_id = set(real_sample_union_id)
    
    # print header first. 
    header = []
    header.append("real_union_id")
    for each in pre_sample_list:
        header.append(each+"_count(pre)")
        header.append(each+"_freq(pre)")
    for each in mock_sample_list:
        header.append(each+"_count(mock)")
        header.append(each+"_freq(mock)")
    for each in real_sample_list:
        header.append(each+"_count(real)")
        header.append(each+"_freq(real)")
    for i in mock_sample_list:
        for j in pre_sample_list:
            header.append("Mock/Pre("+i+"/"+j+")")
    for i in real_sample_list:
        for j in pre_sample_list:
            header.append("Real/Pre(" + i+"/"+j+")")
    for i in real_sample_list:
        for j in mock_sample_list:
            header.append("Real/Mock("+i+"/"+j+")")
    print("{}".format(",".join(header)))

    ##########################
    for each_id in real_sample_union_id:
        print_content = []
        print_content.append(each_id)
        for pre_sample in pre_sample_list:
            print_content.append(big_dict[pre_sample].get(each_id, {'perfect_count':'0'})['perfect_count'])
            print_content.append(big_dict[pre_sample].get(each_id, {'perfect_freq' :'0'})['perfect_freq'])
        for mock_sample in mock_sample_list:
            print_content.append(big_dict[mock_sample].get(each_id, {'perfect_count':'0'})['perfect_count'])
            print_content.append(big_dict[mock_sample].get(each_id, {'perfect_freq' :'0'})['perfect_freq'])
        for real_sample in real_sample_list:
            print_content.append(big_dict[real_sample].get(each_id, {'perfect_count':'0'})['perfect_count'])
            print_content.append(big_dict[real_sample].get(each_id, {'perfect_freq' :'0'})['perfect_freq'])
        # print(print_content)

        # calculate mock/pre:
        foldchange_content = []
        for i in sample_index_dict['mock']:
            for j in sample_index_dict['pre']:
                try:
                    foldchange_content.append(str(float(print_content[i])/float(print_content[j])))
                except ZeroDivisionError:
                    foldchange_content.append('0')
        for i in sample_index_dict['real']:
            for j in sample_index_dict['pre']:
                try:
                    foldchange_content.append(str(float(print_content[i])/float(print_content[j])))
                except ZeroDivisionError:
                    foldchange_content.append('0')
        for i in sample_index_dict['real']:
            for j in sample_index_dict['mock']:
                try:
                    foldchange_content.append(str(float(print_content[i])/float(print_content[j])))
                except ZeroDivisionError:
                    foldchange_content.append('0')
        # print(foldchange_content)
        print("{},{}".format(",".join(print_content), ",".join(foldchange_content)))
        

if __name__ == '__main__':
    main()
