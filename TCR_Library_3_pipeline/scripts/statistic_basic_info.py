""" these modules are used below"""
import sys
import os
import re


def usage():
    """
    docstring for usage
    python statistic_basic_info.py  ./ >basic_info.xls
    """


def addtwodimdict(thedict, key_a, key_b, val):
    ''' this is a function to add two dimetion dict '''
    if key_a in thedict:
        thedict[key_a].update({key_b: val})
    else:
        thedict.update({key_a: {key_b: val}})
    return thedict


def deal_cutadapt_file(file_list):
    """ this function is used for ..."""
    adict = {}
    for afile in file_list:
        file_name = os.path.basename(afile).split(".")[0]
        with open(afile, "r") as f:
            for line in f:
                line = line.rstrip("\n")
                c = re.match(r'(Total reads processed:)\s+(.*)', line)
                if c:
                    # print file_name, c.groups()
                    addtwodimdict(adict, file_name, c.group(1), c.group(2))
                p = re.match(r'(Reads with adapters:)\s+(.*)', line)
                if p:
                    addtwodimdict(adict, file_name, p.group(1), p.group(2))
                    # print file_name, p.groups()
    return adict


def deal_mixcr_log(file_list):
    """ this function is used for ..."""
    adict = {}
    for afile in file_list:
        file_name = os.path.basename(afile).split(".")[0]
        with open(afile, "r") as f:
            for line in f:
                line = line.rstrip("\n")
                a = re.match(r'(Total sequencing reads:)\s+(.*)', line)
                b = re.match(r'(Successfully aligned reads:)\s+(.*)', line)
                c = re.match(r'(TRA chains:)(.*)', line)
                d = re.match(r'(TRB chains:)(.*)', line)
                if a:
                    # print file_name, a.groups()
                    addtwodimdict(adict, file_name, a.group(1), a.group(2))
                if b:
                    addtwodimdict(adict, file_name, b.group(1), b.group(2))
                    # print file_name, b.groups()
                if c:
                    addtwodimdict(adict, file_name, c.group(1), c.group(2))
                    # print file_name, c.groups()
                if d:
                    addtwodimdict(adict, file_name, d.group(1), d.group(2))
                    # print file_name, d.groups()
    return adict


work_dir = sys.argv[1]
result_dir = work_dir + '/' + 'result'
if not os.path.exists(result_dir):
    result_dir = './'

cutadapt_info_list = []
mixcr_report_list = []
for name in os.listdir(result_dir):
    if name.endswith("cutadapt.stats"):
        cutadapt_info_list.append(result_dir + '/' + name)
    if name.endswith(".report"):
        mixcr_report_list.append(result_dir + '/' + name)

if len(cutadapt_info_list) != 0:
    cutadatp_info_dict = deal_cutadapt_file(cutadapt_info_list)

    print("sample\tTotal reads processed\tReads with adapters")
    for k, v in cutadatp_info_dict.items():
        print('{0}\t{1}\t{2}'.format(
            k, v['Total reads processed:'], v['Reads with adapters:']))

if len(mixcr_report_list) != 0:
    mixcr_report_dict = deal_mixcr_log(mixcr_report_list)
    print("""Sample\tTotal sequencing reads\tSuccessfully aligned reads\tTRA chains\tTRB chains""")
    for k, v in mixcr_report_dict.items():
        print('{0}\t{1}\t{2}\t{3}\t{4}'.format(
            k, v.get('Total sequencing reads:',0), v.get('Successfully aligned reads:',0),
            v.get('TRA chains:',0), v.get('TRB chains:',0)))
