import sys
import os
import pandas as pd

files = sys.argv[1:]
for each_file in files:
    filename = os.path.basename(each_file).split(".")[0]
    df = pd.read_table(each_file, header=T)
    print(df)
    pd.melt(df, id_vars="",value_vars=[])

