import sys
import pandas as pd
for each_file in sys.argv[1:]:
    filename = each_file.split()
    data = pd.read_table(each_file, header=T)

