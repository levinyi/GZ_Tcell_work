import sys
import os
import pandas as pd

file_dir = sys.argv[1]
all_file_list = [each for each in os.listdir(file_dir) if each.endswith(".txt")]
#print(all_file_list)

for single_file in all_file_list:
    single_data_frame = pd.read_table(os.path.join(file_dir, single_file), sep=":")
    name = single_file.split(".")[0]
    # change the column names
    single_data_frame.columns = ['Sample', name]

    if single_file == all_file_list[0]:
        all_data_frame = single_data_frame
    else:
        all_data_frame = pd.merge(all_data_frame, single_data_frame, on="Sample")

# print(all_data_frame)
all_data_frame = all_data_frame.T.sort_index()
all_data_frame.to_csv("./summary_report/Summary.prokka.results.xls", index=True,header=False, sep="\t")
