import sys
import pandas as pd

print("Sample.ID\tDetected.Cell.Num\tDetected.TCRpair.Num\tFraction.Reads.in.Cells\tMedian.TRA.UMIs.per.Cell\tMedian.TRB.UMIs.per.Cell")
for each in sys.argv[1:]:
    name = each.split("/")[-3]
    data = pd.read_csv(each)
    for index, row in data.iterrows():
        print("{}\t{}\t{}\t{}\t{}\t{}".format(name, row[0],row[2],row[13],row[14],row[15]))
