import sys,os
import json


def usage():
    print("""
    Usage: 
        python {0}  [file1] [file2]
    
    example: 
        python {0} clonotypes.csv consensus_annotations.csv > xxx/clonotypes_add_consensus_annotations.csv
    
    Update:
    20200529: add RNA-seq data info.
    """.format(os.path.basename(sys.argv[0])))

def read_clonotypes(file1):
	adict = {}
	with open(file1,"r") as f1:
		for line in f1:
			line = line.rstrip("\n")
			if line.startswith("clonotype_id"):
				continue
			clonotype_id,frequency, proportion,cdr3s_aa,cdr3s_nt = line.split(",")
			
			adict.setdefault(clonotype_id,{})['freq'] = frequency
			adict.setdefault(clonotype_id,{})['prop'] = proportion
	return adict


def two_dim_dict(thedict, key_a, key_b, value):
    if key_a in thedict:
        if key_b in thedict[key_a]:
            value = thedict[key_a][key_b] + value
            thedict[key_a].update({key_b: value})
        else:
            thedict[key_a].update({key_b: value})
    else:
        thedict.update({key_a: {key_b: value}})
    return thedict


def read_consensus(file2):
	len_dict = {}
	detail_dict = {}
	with open(file2, "r") as f2:
		for line in f2:
			line = line.rstrip("\n")
			if line.startswith("clonotype_id"):
				continue
			clonotype_id,consensus_id,length,chain,v_gene,d_gene,j_gene,c_gene,full_length,productive,cdr3,cdr3_nt,reads,umis = line.split(",")
			two_dim_dict(len_dict,clonotype_id,chain,1)
			if clonotype_id in detail_dict:
				if chain in detail_dict[clonotype_id]:
					new_chain = chain + "_1"
					detail_dict.setdefault(clonotype_id,{})[new_chain] = [consensus_id,length,chain,v_gene,d_gene,j_gene,c_gene,full_length,productive,cdr3,cdr3_nt,reads,umis]
				else:
					detail_dict.setdefault(clonotype_id,{})[chain] = [consensus_id,length,chain,v_gene,d_gene,j_gene,c_gene,full_length,productive,cdr3,cdr3_nt,reads,umis]
			else:
				detail_dict.setdefault(clonotype_id,{})[chain] = [consensus_id,length,chain,v_gene,d_gene,j_gene,c_gene,full_length,productive,cdr3,cdr3_nt,reads,umis]
	return len_dict, detail_dict


def main():
	# file1 = "GB001002E1L1_TCRseq/clonotypes.csv"
	# file2 = "GB001002E1L1_TCRseq/consensus_annotations.csv"
        if len(sys.argv) != 3:
            usage()
            sys.exit()
        file1 = sys.argv[1]
	file2 = sys.argv[2]
	clonotype_dict = read_clonotypes(file1)
	len_dict, detail_dict = read_consensus(file2)
	# print(json.dumps(detail_dict, indent=4))
	print("clonotype_id,consensus_id,length,chain,v_gene,d_gene,j_gene,c_gene,full_length,productive,cdr3,cdr3_nt,reads,umis,consensus_id,length,chain,v_gene,d_gene,j_gene,c_gene,full_length,productive,cdr3,cdr3_nt,reads,umis,frequency, proportion")
	for each in len_dict:
		if len(len_dict[each]) == 2:
			if len_dict[each]['TRA'] == 1 and len_dict[each]['TRB'] == 1:
				print("{},{},{},{},{}".format(each,",".join(detail_dict[each]['TRA']),",".join(detail_dict[each]['TRB']),clonotype_dict[each]['freq'],clonotype_dict[each]['prop']))
			elif len_dict[each]['TRA'] == 1 and len_dict[each]['TRB'] == 2:
				print("{},{},{},{},{}".format(each+".1",",".join(detail_dict[each]['TRA']),",".join(detail_dict[each]['TRB']),clonotype_dict[each]['freq'],clonotype_dict[each]['prop']))
				print("{},{},{},{},{}".format(each+".2",",".join(detail_dict[each]['TRA']),",".join(detail_dict[each]['TRB_1']),clonotype_dict[each]['freq'],clonotype_dict[each]['prop']))
			elif len_dict[each]['TRA'] == 2 and len_dict[each]['TRB'] == 1:
				print("{},{},{},{},{}".format(each+".1",",".join(detail_dict[each]['TRA']),",".join(detail_dict[each]['TRB']),clonotype_dict[each]['freq'],clonotype_dict[each]['prop']))
				print("{},{},{},{},{}".format(each+".2",",".join(detail_dict[each]['TRA_1']),",".join(detail_dict[each]['TRB']),clonotype_dict[each]['freq'],clonotype_dict[each]['prop']))
			elif len_dict[each]['TRA'] == 2 and len_dict[each]['TRB'] == 2:
				print("{},{},{},{},{}".format(each+".1",",".join(detail_dict[each]['TRA']),",".join(detail_dict[each]['TRB']),clonotype_dict[each]['freq'],clonotype_dict[each]['prop']))
				print("{},{},{},{},{}".format(each+".2",",".join(detail_dict[each]['TRA_1']),",".join(detail_dict[each]['TRB_1']),clonotype_dict[each]['freq'],clonotype_dict[each]['prop']))			
				print("{},{},{},{},{}".format(each+".3",",".join(detail_dict[each]['TRA']),",".join(detail_dict[each]['TRB_1']),clonotype_dict[each]['freq'],clonotype_dict[each]['prop']))
				print("{},{},{},{},{}".format(each+".4",",".join(detail_dict[each]['TRA_1']),",".join(detail_dict[each]['TRB']),clonotype_dict[each]['freq'],clonotype_dict[each]['prop']))			
			# print(each,len_dict[each],len(len_dict[each]))

if __name__ == '__main__':
	main()



