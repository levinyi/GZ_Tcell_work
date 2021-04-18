import os
import sys
import argparse
import numpy as np
import pandas as pd
import configparser

def usage():
    '''
    Usage:
        python3 SnS_Pre_Mock_Real_intergration.py -c SnS_Pre_Mock_Real_intergration.Comp-CR13-CD8T.config  -p Comp-CR13-CD8T -s Pre,Mock,Real

    Updated:
    20200415: updated a new version. add -s parmarter. to calculate average value to a specific group.
    20201222: fix a bug. when config file add space between two samples.
    '''

def _argparse():
    parser = argparse.ArgumentParser(description="This is description")
    parser.add_argument('-c', '--config',  action='store', dest='config', default='config.txt', help="this is config file")
    parser.add_argument('-p', '--project', action='store', dest='project', default='my_project', help='this is project')
    parser.add_argument('-s', '--average', action='store', dest='average', help='default [null].set [Pre,Mock,Real]')
    parser.add_argument('-m', '--min_freq', action='store', dest='min_freq',default='0.0000001', help='min_freq default: 0.0000001')
    return parser.parse_args()


def generate_big_dict(all_sample_list, min_freq,real_sample_list):
    files = os.listdir(path='.')
    big_df = pd.DataFrame()
    new_sample_list = []
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
        df1 = pd.DataFrame(df, columns=['pair','perfect_freq','perfect_count'])
        # replace 0 to min_freq
        df1.loc[df1['perfect_freq']==0,'perfect_freq'] = float(min_freq)
        df1.columns = ['pair', sample+'_freq', sample+'_count']
        if big_df.empty:
            big_df = df1
        else:
            big_df = pd.merge(big_df, df1, on='pair',how='outer')
        new_sample_list.append(sample+'_freq')
    # big_df.to_csv("big_df.raw.xls",sep="\t")
    
    # filter real sample list 
    # real_sample_list2 = [i+"_freq" for i in real_sample_list]
    # big_df = big_df[big_df[real_sample_list2].sum(axis=1) != 0]
    # big_df.to_csv("big_df.filter.xls",sep="\t")
    return big_df, new_sample_list


def check_sample(config_dict):
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
    print("pre_sample : {}\t mock_sample : {}\t Real_sample : {}".format(pre_sample_list, mock_sample_list, real_sample_list))
    return all_sample_list,pre_sample_list,mock_sample_list,real_sample_list


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
    all_sample_list,pre_sample_list,mock_sample_list,real_sample_list = check_sample(config_dict) 
    big_df, new_sample_list = generate_big_dict(all_sample_list, min_freq, real_sample_list)
    #print(big_df)

    if parser.average:
        ave_set = parser.average.rstrip(",").split(",")  # [Pre,Mock,Real]
        print("since you set average group is {}".format(ave_set))
        for each in ave_set: # [Pre,Mock,Real]
            if each == 'Pre':
                pre_sample_incolumn = [i+"_freq" for i in pre_sample_list]
                big_df['pre_ave_freq'] = big_df[pre_sample_incolumn].mean(axis=1)
                pre_sample_list = ["pre_ave"]
            elif each == 'Mock':
                mock_sample_incolumn = [i+"_freq" for i in mock_sample_list]
                big_df['mock_ave_freq'] = big_df[mock_sample_incolumn].mean(axis=1)
                mock_sample_list = ["mock_ave"]
            elif each == 'Real':
                real_sample_incolumn = [i+"freq" for i in mock_sample_list]
                big_df['real_ave_freq'] = big_df[real_sample_incolumn].mean(axis=1)
                real_sample_list = ["real_ave"]
        # print(big_df)
        # big_df.to_csv("big_df.cal_ave.csv",sep="\t")
    print(pre_sample_list,mock_sample_list,real_sample_list)

    ##############################################
    for pre_sample in pre_sample_list:
        for mock_sample in mock_sample_list:
            big_df['Mock/Pre({}/{})'.format(mock_sample,pre_sample)] = big_df[mock_sample+"_freq"]/big_df[pre_sample+"_freq"]
        for real_sample in real_sample_list:
            big_df['Real/Pre({}/{})'.format(real_sample,pre_sample)] = big_df[real_sample+"_freq"]/big_df[pre_sample+"_freq"]
    for moc_sample in mock_sample_list:
        for real_sample in real_sample_list:
            big_df['Real/Mock({}/{})'.format(real_sample,mock_sample)] = big_df[real_sample+"_freq"]/big_df[mock_sample+"_freq"]
    #print(big_df)
    ####################################################################
    # edit header for. 
    big_df = big_df.rename(columns={'pair':'TCR_id(Real_samples_union)'})
    #print(big_df)
    big_df.to_csv(parser.project+".csv",sep=",",index=False)
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
