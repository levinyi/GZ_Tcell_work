import sys
import pandas as pd
vdj_db = "/cygene/software/biosoftware/vdjdb-db/database/vdjdb.txt"
vdj_db_table = pd.read_table(vdj_db, sep="\t",).fillna(value="NA")

my_data = sys.argv[1]
my_data = "consensus_annotations.csv"
my_data_table = pd.read_table(my_data, sep=",").fillna(value="NA")

new  = my_data_table.loc[my_data_table['cdr3'].isin(vdj_db_table['cdr3'].tolist()),'cdr3']
