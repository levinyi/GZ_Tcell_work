import sys
import os
import pandas as pd

file_dir = sys.argv[1]
all_file_dir = [each for each in os.listdir(os.path.abspath(file_dir)) if os.path.isdir(os.path.join(os.path.abspath(file_dir),each))]
all_file_list =[each+'/'+each +'.megahit-quast-report/report.tsv' for each in all_file_dir]

for single_file in all_file_list:
    single_data_frame = pd.read_table(os.path.join(file_dir,single_file), sep="\t")
    name = single_file.split("/")[0]
    # change the column names
    single_data_frame.columns = ['Sample', name]
    if single_file == all_file_list[0]:
        all_data_frame = single_data_frame
    else:
        all_data_frame = pd.merge(all_data_frame, single_data_frame, on="Sample")

all_data_frame = all_data_frame.T.sort_index()
all_data_frame.to_csv("./summary_report/Summary.megahit.quest.results.xls", index=True, header=False, sep="\t")

