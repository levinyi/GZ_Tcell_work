import sys
import os

def addtwodimdict(thedict, key_a, key_b, val):
    ''' this is a function to add two dimetion dict '''
    if key_a in thedict:
        thedict[key_a].update({key_b: val})
    else:
        thedict.update({key_a: {key_b: val}})
    return thedict


cell_number = sys.argv[1]

cell_number_dict = {}
with open(cell_number, "r") as f:
    for line in f:
        line = line.rstrip("\n")
        sample, cell = line.split()
        cell_number_dict[sample] = cell

sample_list = []
folder_list = []
sample_umi_dict = {}
sample_cell_dict = {}
sample_mole_dict = {}

folder = sys.argv[2]
for each in os.listdir(folder):
    if os.path.isdir(each):
        #print("yes {}".format(each))
        folder_list.append(each)
        sub_dir = os.path.abspath(folder) + '/' + each
        for each_file in os.listdir(sub_dir):
            if 'pair.count' in each_file:
                if os.path.isfile(sub_dir + '/' + each_file):
                    sample_name = each_file.split(".")[0]
                    sample_list.append(sample_name)
                    with open(sub_dir + '/' + each_file, "r") as f:
                        for line in f:
                            line = line.rstrip("\n")
                            if line.startswith("total_reads"):
                                continue
                            c = line.split() # total_reads match_reads perfect_reads percent_perfect_reads umi_total umi_perfect percent_perfect_umi
                            umi_total = c[4] 
                            cell_total = float(cell_number_dict[sample_name])
                            molecule_per_cell = float(umi_total) / cell_total
                            #print("{}\t{}\t{}\t{}\t{}".format(each, sample_name, umi_total, cell_total, molecule_per_cell))
                            addtwodimdict(sample_umi_dict, sample_name, each, umi_total)
                            addtwodimdict(sample_mole_dict, sample_name, each, molecule_per_cell)
a = []
for i in folder_list:
    a.append("umi_total_in_"+i)
    a.append("molecule_per_cell_in_"+i)
print("sample,cell_number,{}".format(",".join(a)))
unique_sample_list = cell_number_dict.keys()
for sample in unique_sample_list:
    printed_content = []
    printed_content.append(sample)
    printed_content.append(cell_number_dict[sample])
    for folder in folder_list:
        printed_content.append(sample_umi_dict[sample][folder])
        printed_content.append(sample_mole_dict[sample][folder])
    print("{}".format(",".join(map(str,printed_content))))
