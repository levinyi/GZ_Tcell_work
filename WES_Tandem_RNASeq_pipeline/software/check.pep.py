import sys
import re
from Bio import SeqIO
from multiprocessing import Pool


def deal_fasta(fastfile):
    fasta_dict = {}
    for record in SeqIO.parse(fastfile, "fasta"):
        if str(record.id).endswith("_1"):
            record_id = str(record.id)[:-2]
        else:
            record_id = str(record.id).split(".")[0]
        #match = re.match(r'(\w+_\d+)\_?1?', record_id)
        #fasta_dict[match.group(1)] = str(record.seq).upper()
        fasta_dict[record_id] = str(record.seq)
    return fasta_dict


def addtwodimdict(thedict, key_a, key_b, val):
    ''' this is a function to add two dimetion dict '''
    if key_a in thedict:
        thedict[key_a].update({key_b: val})
    else:
        thedict.update({key_a: {key_b: val}})
    return thedict


def find_mutation_pos(wildtype_seq, mutation_seq):
    """docstring for find_mutation"""
    if len(wildtype_seq) < len(mutation_seq):
        for i in range(0, len(wildtype_seq)):
            mutation_pos = i+1
            wildtype_aa = wildtype_seq[i]
            mutation_aa = mutation_seq[i]
            if wildtype_aa != mutation_aa:
                return wildtype_aa, mutation_pos, mutation_aa
        return "*", mutation_pos, mutation_seq[mutation_pos]
    else: # stop loss 
        for i in range(0, len(mutation_seq)):
            mutation_pos = i+1
            wildtype_aa = wildtype_seq[i]
            mutation_aa = mutation_seq[i]
            if wildtype_aa != mutation_aa:
                return wildtype_aa, mutation_pos, mutation_aa
    return wildtype_seq[mutation_pos], mutation_pos, "Q"


def find_mutation_pos_list(wildtype_seq, mutation_seq):
    alist = []
    if len(wildtype_seq) <= len(mutation_seq):
        range_length = len(wildtype_seq)
    else:
        range_length = len(mutation_seq)

    for i in range(0, range_length):
        if wildtype_seq[i] != mutation_seq[i]:
            mutation_pos = i+1
            wildtype_aa = wildtype_seq[i]
            mutation_aa = mutation_seq[i]
            alist.append([mutation_pos, wildtype_aa, mutation_aa])
    return alist


def deal_STF(stf_file):
    trans_id_dict = {}
    with open(stf_file, "r") as f:
        for line in f:
            line = line.rstrip("\n")
            trans_id_dict.setdefault(line.split("\t")[7], []).append(line)
    return trans_id_dict


def main():
    '''
    with open(my_data, "r") as f:
        for line in f:
            line = line.rstrip("\n")
            if line.startswith("OM_ID"):
                continue
            if len(line.split("\t")) == 11:
                # STF file
                Chromosome, Position_start, Position_end,Ref,Alt, annotation,\
                Gene,Transcript_Ref,exon,Coding_DNA_change,AA_change = line.split("\t")
            else:
                OM_ID, Gene, Transcript_Ref, Classification, Alternation_type1,\
                    Alternation_type2, Chromosome, Position_start, Position_end,\
                    Genomic_DNA_change, Coding_DNA_change, AA_change, \
                    Vaf, Mean_coverage = line.split("\t")
    '''
    p = Pool(2)
    # hg19pep = "hg19.pep.fa"
    hg19pep = sys.argv[1]
    my_data = sys.argv[2]
    mutation_pep = sys.argv[3]
    half_minigene_length = 14
    # mutation_pep = "genome_replaced_mutation.pep.fa"
    # my_data = "data1.xls"
    p1 = p.apply_async(deal_fasta, args=(hg19pep,))
    p2 = p.apply_async(deal_fasta, args=(mutation_pep,))
    p.close()
    p.join()
    hg19_fasta_dict = p1.get()
    mutation_fasta_dict = p2.get()
    trans_id_dict = deal_STF(my_data)
    # print(len(hg19_fasta_dict))
    # print(len(mutation_fasta_dict))

    for trans_id in trans_id_dict: 
        if trans_id in hg19_fasta_dict and trans_id in mutation_fasta_dict:
            wildtype_seq = hg19_fasta_dict[trans_id]
            mutation_seq = mutation_fasta_dict[trans_id]
        else:
            print("{}\tNA\tNA\tNA\tNA".format("\t".join(trans_id_dict[trans_id])))
            continue

        if len(trans_id_dict[trans_id]) > 1:
            print("this trans id has more than 2 mutation sites ", trans_id)
            trans_mut_info_list = trans_id_dict[trans_id]
            mutation_infolist = find_mutation_pos_list(wildtype_seq, mutation_seq)  # ['line1', 'line2']
            sorted_mutation_list = sorted(mutation_infolist)
            #sorted_mutation_list2 = ["p.{}{}{}".format(each[1],each[0],each[2]) for each in sorted_mutation_list]
            
            new_dict = {each.split("\t")[10]: each for each in trans_mut_info_list}
            new_list = new_dict.keys()
            sorted_trans_mut_infolist = sorted(new_list,reverse=True)
            '''
            if len(sorted_mutation_list) > len(sorted_trans_mut_infolist):
                range_length = len(sorted_trans_mut_infolist)
            else:
                range_length = len(sorted_mutation_list)
            '''
            print("I found these mutation aa change:", sorted_mutation_list, len(sorted_mutation_list))
            print("you reported these mutation aa change", sorted_trans_mut_infolist, len(sorted_trans_mut_infolist))
            for i in range(0, len(sorted_trans_mut_infolist)):
                try:
                    mutation_pos, wildtype_aa, mutation_aa = sorted_mutation_list[i][0],sorted_mutation_list[i][1],sorted_mutation_list[i][2]
                except:
                    print("{}\tNA\tNO\tNA\tNA".format(trans_id_dict[trans_id][i]))
                    continue
                # print(type(mutation_pos), mutation_pos, wildtype_aa, mutation_aa)
                find_aa_change = "p.{}{}{}".format(wildtype_aa, mutation_pos, mutation_aa)
                reported_aa_change = sorted_trans_mut_infolist[i]
                if find_aa_change == reported_aa_change:   # mutation_pos, wildtype_aa, mutation_aa
                    print("yes")
                    if 0 < mutation_pos <= half_minigene_length:
                        print("yes1")
                        print("{}\t{}\tYES\t#{}\t#{}".format(
                            trans_id_dict[trans_id][i], 
                            find_aa_change,
                            wildtype_seq[:mutation_pos + half_minigene_length], 
                            mutation_seq[:mutation_pos + half_minigene_length],
                            )
                        )
                    elif len(mutation_seq) >= mutation_pos > len(mutation_seq) - half_minigene_length:
                        print("yes2")

                        print("{}\t{}\tYES\t{}#\t{}#".format(
                            trans_id_dict[trans_id][i], 
                            find_aa_change, 
                            wildtype_seq[mutation_pos -half_minigene_length-1:], 
                            mutation_seq[mutation_pos -half_minigene_length-1:],
                            )
                        )
                    else:
                        print("yes3")

                        print("{}\t{}\tYES\t{}\t{}".format(
                            trans_id_dict[trans_id][i], 
                            find_aa_change, 
                            wildtype_seq[mutation_pos-half_minigene_length-1:mutation_pos+half_minigene_length], 
                            mutation_seq[mutation_pos-half_minigene_length-1:mutation_pos+half_minigene_length],
                            )
                        )
                else:
                    #print("no")
                    if 0 < mutation_pos <= half_minigene_length:
                        print("no1")
                        print("{}\t{}\tNO\t#{}\t#{}".format(
                            trans_id_dict[trans_id][i], 
                            find_aa_change, 
                            wildtype_seq[:mutation_pos+half_minigene_length], 
                            mutation_seq[:mutation_pos+half_minigene_length],
                            )
                        )
                    elif len(mutation_seq) > mutation_pos > len(mutation_seq) - half_minigene_length:
                        print("no2")
                        print("{}\t{}\tNO\t{}#\t{}#".format(
                            trans_id_dict[trans_id][i], 
                            find_aa_change, 
                            wildtype_seq[mutation_pos -half_minigene_length-1:], 
                            mutation_seq[mutation_pos -half_minigene_length-1:],
                            )
                        )
                    else:
                        print("no3")
                        print("{}\t{}\tNO\t{}\t{}".format(
                            trans_id_dict[trans_id][i], 
                            find_aa_change,
                            wildtype_seq[mutation_pos-half_minigene_length-1:mutation_pos+half_minigene_length], 
                            mutation_seq[mutation_pos-half_minigene_length-1:mutation_pos+half_minigene_length],
                            )
                        )
        else:
            if wildtype_seq == mutation_seq:
                print("{}\t-\tsynonymous\t-\t-".format(line))
            else:
                line = trans_id_dict[trans_id][0]
                # print(line)
                AA_change = line.split("\t")[10]
                wildtype_aa, mutation_pos, mutation_aa = find_mutation_pos(wildtype_seq, mutation_seq)
                if "p.{}{}{}".format(wildtype_aa, mutation_pos, mutation_aa) == AA_change:
                    if 0 < mutation_pos <= half_minigene_length:
                        print("{}\tp.{}{}{}\tYES\t#{}\t#{}".format(
                            line, 
                            wildtype_aa, mutation_pos, mutation_aa, 
                            wildtype_seq[:mutation_pos + half_minigene_length], 
                            mutation_seq[:mutation_pos + half_minigene_length],
                            )
                        )
                    elif len(mutation_seq) >= mutation_pos > len(mutation_seq) - half_minigene_length:
                        print("{}\tp.{}{}{}\tYES\t{}#\t{}#".format(
                            line, 
                            wildtype_aa, mutation_pos, mutation_aa, 
                            wildtype_seq[mutation_pos -half_minigene_length-1:], 
                            mutation_seq[mutation_pos -half_minigene_length-1:],
                            )
                        )
                    else:
                        print("{}\tp.{}{}{}\tYES\t{}\t{}".format(
                            line, 
                            wildtype_aa, mutation_pos, mutation_aa, 
                            wildtype_seq[mutation_pos-half_minigene_length-1:mutation_pos+half_minigene_length], 
                            mutation_seq[mutation_pos-half_minigene_length-1:mutation_pos+half_minigene_length],
                            )
                        )
                else:
                    if half_minigene_length > mutation_pos > 0:
                        print("{}\tp.{}{}{}\tNO\t#{}\t#{}".format(
                            line, 
                            wildtype_aa, mutation_pos, mutation_aa, 
                            wildtype_seq[:mutation_pos+half_minigene_length], 
                            mutation_seq[:mutation_pos+half_minigene_length],
                            )
                        )
                    elif len(mutation_seq) > mutation_pos > len(mutation_seq) - half_minigene_length:
                        print("{}\tp.{}{}{}\tNO\t{}#\t{}#".format(
                            line, 
                            wildtype_aa, mutation_pos, mutation_aa, 
                            wildtype_seq[mutation_pos -half_minigene_length-1:], 
                            mutation_seq[mutation_pos -half_minigene_length-1:],
                            )
                        )
                    else:
                        print("{}\tp.{}{}{}\tNO\t{}\t{}".format(
                            line, 
                            wildtype_aa, mutation_pos, mutation_aa,
                            wildtype_seq[mutation_pos-half_minigene_length-1:mutation_pos+half_minigene_length], 
                            mutation_seq[mutation_pos-half_minigene_length-1:mutation_pos+half_minigene_length],
                            )
                        )
                    # elif reported_mutation_aa in diff_pos_dict[reported_wildtype_aa].values():


if __name__ == '__main__':
    main()
