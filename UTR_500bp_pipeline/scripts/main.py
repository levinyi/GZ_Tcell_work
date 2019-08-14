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

def main():
    parser = _argparse()
    cf = ConfigParser.ConfigParser(allow_no_value=True)
    cf.read(parser.config)

    # check config file
    if not cf.has_section('Config'):
        os.exit("your config file is not correct.")

    # read config file
    adapter1 = cf.get('Config', 'adapter1')
    adapter2 = cf.get('Config', 'adapter2')
    database = cf.get('Config', 'database')
    scripts = cf.get('Config', 'scripts')
    project_name = cf.get("Config", "project_name")
    merge_analysis = cf.get("Config", "merge_analysis")
    library_type = cf.get("Config", "library_type")

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
    shell_name = shell_dir + '/work.sh'
    step1_shell = open(shell_name, "w")
    rawdata_dict = deal_rawdata(parser.data, data_dir)
