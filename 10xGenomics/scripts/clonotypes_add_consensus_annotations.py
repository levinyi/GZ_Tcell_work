import sys,os
import json


def usage():
    print("""
    Usage: 
        python {0}  [file1] [file2]
    
    example: 
        python {0} clonotypes.csv consensus_annotations.csv > xxx/clonotypes_add_consensus_annotations.csv
        python3 {0} > clonotypes_add_consensus_annotations.csv
    
    Update:
    20200709: Add an example for no input file situation.
    20200605: optimized code: if there is no input file, default search  clonotypes.csv and consensus_annotations.csv in current directory.
    20200529: add RNA-seq data info.
    """.format(os.path.basename(sys.argv[0])))


def read_clonotypes(file1):
    adict = {}
    with open(file1,"r") as f1:
        for line in f1:
            line = line.rstrip("\n")
            if line.startswith("clonotype_id"):
                continue
            clonotype_id,frequency, proportion,cdr3s_aa,cdr3s_nt,inkt_evidence,mait_evidence = line.split(",")
            # for old clonotype file
            # clonotype_id,frequency,proportion,cdr3s_aa,cdr3s_nt = line.split(",")
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
            a = line.split(",")
            clonotype_id = a[0]
            chain = a[3]

            two_dim_dict(len_dict,clonotype_id,chain,1)
            if clonotype_id in detail_dict:
                if chain in detail_dict[clonotype_id]:
                    new_chain = chain + "_1"
                    detail_dict.setdefault(clonotype_id,{})[new_chain] = a[1:]
                else:
                    detail_dict.setdefault(clonotype_id,{})[chain] = a[1:]
            else:
                detail_dict.setdefault(clonotype_id,{})[chain] = a[1:]
    return len_dict, detail_dict


def main():
	# file1 = "GB001002E1L1_TCRseq/clonotypes.csv"
	# file2 = "GB001002E1L1_TCRseq/consensus_annotations.csv"
        if len(sys.argv) == 1:
            file1 = "./clonotypes.csv"
            file2 = "./consensus_annotations.csv"
        elif len(sys.argv) == 3:
            file1 = sys.argv[1]
            file2 = sys.argv[2]
        else:
            usage()
            sys.exit()
        output = open("clonotypes_add_consensus_annotations.csv", "w")
        
        clonotype_dict = read_clonotypes(file1)
        len_dict, detail_dict = read_consensus(file2)
	# print(json.dumps(detail_dict, indent=4))
        output.write("clonotype_id,consensus_id,length,chain,v_gene,d_gene,j_gene,c_gene,full_length,productive,fwr1,fwr1_nt,cdr1,cdr1_nt,fwr2,fwr2_nt,cdr2,cdr2_nt,fwr3,fwr3_nt,cdr3,cdr3_nt,fwr4,fwr4_nt,reads,umis,v_start,v_end,v_end_ref,j_start,j_start_ref,j_end,fwr1_start,fwr1_end,cdr1_start,cdr1_end,fwr2_start,fwr2_end,cdr2_start,cdr2_end,fwr3_start,fwr3_end,cdr3_start,cdr3_end,fwr4_start,fwr4_end,consensus_id,length,chain,v_gene,d_gene,j_gene,c_gene,full_length,productive,fwr1,fwr1_nt,cdr1,cdr1_nt,fwr2,fwr2_nt,cdr2,cdr2_nt,fwr3,fwr3_nt,cdr3,cdr3_nt,fwr4,fwr4_nt,reads,umis,v_start,v_end,v_end_ref,j_start,j_start_ref,j_end,fwr1_start,fwr1_end,cdr1_start,cdr1_end,fwr2_start,fwr2_end,cdr2_start,cdr2_end,fwr3_start,fwr3_end,cdr3_start,cdr3_end,fwr4_start,fwr4_end,frequency,proportion\n")
        # output.write("clonotype_id,consensus_id,length,chain,v_gene,d_gene,j_gene,c_gene,full_length,productive,cdr3,cdr3_nt,reads,umis,consensus_id,length,chain,v_gene,d_gene,j_gene,c_gene,full_length,productive,cdr3,cdr3_nt,reads,umis,frequency,proportion\n")
        for each in len_dict:
            if len(len_dict[each]) == 2:
                if len_dict[each]['TRA'] == 1 and len_dict[each]['TRB'] == 1:
                    output.write("{},{},{},{},{}\n".format(each,",".join(detail_dict[each]['TRA']),",".join(detail_dict[each]['TRB']),clonotype_dict[each]['freq'],clonotype_dict[each]['prop']))
                elif len_dict[each]['TRA'] == 1 and len_dict[each]['TRB'] == 2:
                    output.write("{},{},{},{},{}\n".format(each+".1",",".join(detail_dict[each]['TRA']),",".join(detail_dict[each]['TRB']),clonotype_dict[each]['freq'],clonotype_dict[each]['prop']))
                    output.write("{},{},{},{},{}\n".format(each+".2",",".join(detail_dict[each]['TRA']),",".join(detail_dict[each]['TRB_1']),clonotype_dict[each]['freq'],clonotype_dict[each]['prop']))
                elif len_dict[each]['TRA'] == 2 and len_dict[each]['TRB'] == 1:
                    output.write("{},{},{},{},{}\n".format(each+".1",",".join(detail_dict[each]['TRA']),",".join(detail_dict[each]['TRB']),clonotype_dict[each]['freq'],clonotype_dict[each]['prop']))
                    output.write("{},{},{},{},{}\n".format(each+".2",",".join(detail_dict[each]['TRA_1']),",".join(detail_dict[each]['TRB']),clonotype_dict[each]['freq'],clonotype_dict[each]['prop']))
                elif len_dict[each]['TRA'] == 2 and len_dict[each]['TRB'] == 2:
                    output.write("{},{},{},{},{}\n".format(each+".1",",".join(detail_dict[each]['TRA']),",".join(detail_dict[each]['TRB']),clonotype_dict[each]['freq'],clonotype_dict[each]['prop']))
                    output.write("{},{},{},{},{}\n".format(each+".2",",".join(detail_dict[each]['TRA_1']),",".join(detail_dict[each]['TRB_1']),clonotype_dict[each]['freq'],clonotype_dict[each]['prop']))			
                    output.write("{},{},{},{},{}\n".format(each+".3",",".join(detail_dict[each]['TRA']),",".join(detail_dict[each]['TRB_1']),clonotype_dict[each]['freq'],clonotype_dict[each]['prop']))
                    output.write("{},{},{},{},{}\n".format(each+".4",",".join(detail_dict[each]['TRA_1']),",".join(detail_dict[each]['TRB']),clonotype_dict[each]['freq'],clonotype_dict[each]['prop']))			
                    # print(each,len_dict[each],len(len_dict[each]))
        output.close()

if __name__ == '__main__':
	main()

