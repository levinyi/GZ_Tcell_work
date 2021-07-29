import os
import sys
import argparse
import numpy as np
import pandas as pd
import configparser

def usage():
    '''
    Usage:
        python3 SnS_Pre_Mock_Real_intergration.v3.py -c SnS_Pre_Mock_Real_intergration.Comp-CR13-CD8T.config  -p Comp-CR13-CD8T -s Pre,Mock,Real

    Updated:
    20210729: fix a bug when mock is empty in config file. will report KeyError.
    20200419: fix some bugs.
    20200415: updated a new version. add -s parmarter. to calculate average value to a specific group.
    20201222: fix a bug. when config file add space between two samples.

    Author: shiyi@rootpathgx.com
    '''

def _argparse():
    parser = argparse.ArgumentParser(description="This is description")
    parser.add_argument('-c', '--config',  action='store', dest='config', default='config.txt', help="this is config file")
    parser.add_argument('-p', '--project', action='store', dest='project', default='my_project', help='this is project')
    parser.add_argument('-s', '--average', action='store', dest='average', help='default [null].set [Pre,Mock,Real]')
    parser.add_argument('-m', '--min_freq', action='store', dest='min_freq',default='0.0000001', help='min_freq default: 0.0000001')
    return parser.parse_args()


def generate_big_dict(all_sample_list, min_freq,real_sample_list,sample_dict):
    # input file sample_dict{sample:pre,sample1:mock,sample3:real}
    # print(sample_dict)
    files = os.listdir(path='.')
    big_df = pd.DataFrame()
    new_sample_list = []
    new_sample_dict = {}
    pre_sample_list = []
    mock_sample_list = []
    real_sample_list = []
    for sample in all_sample_list:
        sample_file = [each for each in files if each.startswith(sample+'.pair.acc_freq')]
        print("detected file: {}".format(sample_file))
        if len(sample_file) > 1:
            sys.exit("Error! Sample file duplicated: {}".format(sample_file))
        try:
            sample_file = sample_file[0]
        except IndexError:
            sys.exit("Error: Wrong sample name in your config file: {}. Could not find {}.pair.acc_freq*.txt".format(sample,sample))
        df = pd.read_table(sample_file, sep="\t").fillna(value=min_freq)
        # print(df.head())
        df1 = pd.DataFrame(df, columns=['pair','perfect_count','perfect_freq'])
        # print(df1)
        # replace 0 to min_freq in freq columns.
        df1['perfect_freq'] = df1['perfect_freq'].replace(0,float(min_freq))
        # replace column name : add sample name to each column.
        df1.columns = ['pair', sample+'_count({})'.format(sample_dict[sample]), sample+'_freq({})'.format(sample_dict[sample])]
        new_sample_list.append(sample+'_freq({})'.format(sample_dict[sample]))
        new_sample_dict.setdefault(sample_dict[sample],[]).append(sample+'_freq({})'.format(sample_dict[sample]))
        if sample_dict[sample] == "Pre":
            pre_sample_list.append(sample+'_freq({})'.format(sample_dict[sample]))
        if sample_dict[sample] == "Mock":
            mock_sample_list.append(sample+'_freq({})'.format(sample_dict[sample]))
        if sample_dict[sample] == "Real":
            real_sample_list.append(sample+'_freq({})'.format(sample_dict[sample]))
        # print(df1.head())
        if big_df.empty:
            big_df = df1
        else:
            big_df = pd.merge(big_df, df1, on='pair',how='outer')
    # print(big_df)
    # replace NAN to min_freq in freq columns 
    big_df[new_sample_list] = big_df[new_sample_list].fillna(value=float(min_freq))
    # replace NAN to 0 in count columns.
    big_df = big_df.fillna(0)
    # print(big_df)
    print("new sample list: {}".format(new_sample_list))
    # big_df.to_csv("big_df.raw.xls",sep="\t")
    # real_sample_list2 = [i+"_freq" for i in real_sample_list]
    # big_df = big_df[big_df[real_sample_list2].sum(axis=1) != 0]
    # big_df.to_csv("big_df.filter.xls",sep="\t")
    # print(big_df.head())
    return big_df, new_sample_list, new_sample_dict,pre_sample_list, mock_sample_list, real_sample_list


def check_sample(config_dict):
    sample_dict = {}
    if len(config_dict['pre_samples']) == 0:
        pre_sample_list = []
    else:
        pre_sample_list = config_dict['pre_samples'].replace(" ","").rstrip(",").split(",")
        for each_sample in pre_sample_list:
            sample_dict[each_sample] = 'Pre'
    if len(config_dict['mock_samples']) == 0:
        mock_sample_list = []
    else:
        mock_sample_list = config_dict['mock_samples'].replace(" ","").rstrip(",").split(",")
        for each_sample in mock_sample_list:
            sample_dict[each_sample] = 'Mock'
    if len(config_dict['real_samples']) == 0:
        real_sample_list = []
    else:
        real_sample_list = config_dict['real_samples'].replace(" ","").rstrip(",").split(",")
        for each_sample in real_sample_list:
            sample_dict[each_sample] = 'Real'
    all_sample_list = pre_sample_list + mock_sample_list + real_sample_list
    print("your input sample names are: {}".format(all_sample_list))
    print("pre_sample : {}\t mock_sample : {}\t Real_sample : {}".format(pre_sample_list, mock_sample_list, real_sample_list))
    return all_sample_list,pre_sample_list,mock_sample_list,real_sample_list,sample_dict


def main():
    parser = _argparse()
    min_freq = parser.min_freq

    cf = configparser.ConfigParser(allow_no_value=False, )
    cf.read(parser.config)
    if not cf.has_section('data'):
        os.exit("Error: your config file is not correct.")

    config_dict = {
        'pre_samples' : cf.get('data','Pre'),
        'mock_samples' : cf.get('data', 'Mock'),
        'real_samples' : cf.get('data', 'Real'),
        }
    all_sample_list,pre_sample_list,mock_sample_list,real_sample_list, sample_dict = check_sample(config_dict) 
    big_df, new_sample_list, new_sample_dict ,pre_sample_list, mock_sample_list, real_sample_list = generate_big_dict(all_sample_list, min_freq, real_sample_list,sample_dict)
    # print(big_df)

    # print("new_sample_dict: {}".format(new_sample_dict))
    ##############################################
    if parser.average:
        ave_set = parser.average.rstrip(",").split(",")  # [Pre,Mock,Real]
        print("since you set average group is {}".format(ave_set))
        if 'Pre' in ave_set:  # [Pre, Mock, Real]
            big_df['pre-ave_freq'] = big_df[new_sample_dict['Pre']].mean(axis=1)
            pre_sample_list = ["pre-ave_freq"]
        else:
            pre_sample_list = new_sample_dict.get('Pre',[])
        if 'Mock' in ave_set:
            big_df['mock-ave_freq'] = big_df[new_sample_dict['Mock']].mean(axis=1)
            mock_sample_list = ["mock-ave_freq"]
        else:
            mock_sample_list = new_sample_dict.get('Mock',[])
        if 'Real' in ave_set:
            big_df['real-ave_freq'] = big_df[new_sample_dict['Real']].mean(axis=1)
            real_sample_list = ["real-ave_freq"]
        else:
            real_sample_list = new_sample_dict.get('Real',[])
        # print(big_df.head())
    print("pre_samplelist,mocksamplilist,real_samplelist {} {} {} ".format(pre_sample_list,mock_sample_list,real_sample_list))

    ##############################################
    # calculate mock/pre real/pre  real/mock
    for pre_sample in pre_sample_list:
        for mock_sample in mock_sample_list:
            big_df['Mock/Pre({}/{})'.format(mock_sample.split("_")[0],pre_sample.split("_")[0])] = big_df[mock_sample]/big_df[pre_sample]
    for pre_sample in pre_sample_list:
        for real_sample in real_sample_list:
            big_df['Real/Pre({}/{})'.format(real_sample.split("_")[0],pre_sample.split("_")[0])] = big_df[real_sample]/big_df[pre_sample]
    for mock_sample in mock_sample_list:
        for real_sample in real_sample_list:
            big_df['Real/Mock({}/{})'.format(real_sample.split("_")[0],mock_sample.split("_")[0])] = big_df[real_sample]/big_df[mock_sample]
    # print(big_df)
    ####################################################################
    # edit header for. 
    big_df = big_df.rename(columns={'pair':'TCR_id(Real_samples_union)'})
    #print(big_df)
    big_df.to_csv(parser.project+".csv",sep=",",index=False)
    print("all done!")


if __name__ == '__main__':
    main()
