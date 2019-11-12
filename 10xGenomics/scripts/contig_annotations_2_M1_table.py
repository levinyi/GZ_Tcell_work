import sys



def addtwodimdict(thedict, key_a, key_b, val):
    ''' this is a function to add two dimetion dict '''
    if key_a in thedict:
        thedict[key_a].update({key_b: val})
    else:
        thedict.update({key_a:{key_b: val}})
    return thedict

def deal_file(afile):
    TCR_dict = {}
    with open(afile, "r") as f:
        for line in f:
            line = line.rstrip("\n")
            # clonotype_id,consensus_id,length,chain,v_gene,d_gene,j_gene,c_gene,full_length,productive,cdr3,cdr3_nt,reads,umis = line.split(",")
            barcode,is_cell,contig_id,high_confidence,length,chain,v_gene,d_gene,j_gene,c_gene,full_length,productive,cdr3,cdr3_nt,reads,umis,raw_clonotype_id,raw_consensus_id = line.split(",")
            if is_cell == "True" and cdr3 != "None" and productive == "True":
                if chain =="TRA":
                    # clonotype = v_gene.replace("/","")+cdr3+j_gene
                    addtwodimdict(TCR_dict, barcode, 'TRAJ', j_gene)
                    addtwodimdict(TCR_dict, barcode, 'CDR3Alpha', cdr3)
                    addtwodimdict(TCR_dict, barcode, 'TRAV', v_gene.replace("/",""))
                    addtwodimdict(TCR_dict, barcode, 'Clone_A_reads', reads)
                    addtwodimdict(TCR_dict, barcode, 'Clone_A_umi', umis)
                elif chain == "TRB":
                    # clonotype = v_gene.replace("/","")+cdr3+j_gene
                    addtwodimdict(TCR_dict, barcode, 'TRBV',v_gene.replace("/",""))
                    addtwodimdict(TCR_dict, barcode, 'CDR3Beta',cdr3)
                    addtwodimdict(TCR_dict, barcode, 'TRBJ',j_gene)
                    addtwodimdict(TCR_dict, barcode, 'Clone_B_reads', reads)
                    addtwodimdict(TCR_dict, barcode, 'Clone_B_umi', umis)
                else:
                    print("other vgene")
    return TCR_dict

if __name__ == '__main__':
    contig_file = sys.argv[1]
    TCR_dict = deal_file(contig_file)
    for b in TCR_dict:
        print b,TCR_dict[b].get('TRAV',None),TCR_dict[b].get('CDR3Alpha',None),TCR_dict[b].get('TRAJ',None),\
        TCR_dict[b].get('TRBV',None),TCR_dict[b].get('CDR3Beta',None),TCR_dict[b].get('TRBJ',None),\
        TCR_dict[b].get('Clone_A_reads',None),TCR_dict[b].get('Clone_A_umi',None),TCR_dict[b].get('Clone_B_reads',None),TCR_dict[b].get('Clone_B_umi',None)
