import sys
from Bio import SeqIO
import re
import argparse


def _argparse():
    parser = argparse.ArgumentParser(description="This is description")
    parser.add_argument('-c', '--cds', action='store', dest='cds_file', help="this is config file")
    parser.add_argument('-g', '--gtf', action='store', dest='gtf_file', help="gtf file used for find (+/-) chains.")
    parser.add_argument('-l', '--indel', action='store', dest='indel_file', help="this is data list file")
    parser.add_argument('-o', '--output', action='store', dest='output_file', help="this is data list file")
    return parser.parse_args()


def deal_fasta(fastfile):
    fasta_dict = {}
    for record in SeqIO.parse(fastfile, "fasta"):
        match = re.match(r'(\w+_\d+)\.\d+_?1?', str(record.id))
        fasta_dict[match.group(1)] = str(record.seq)
    return fasta_dict


def deal_gtf(gtf_file):
    """ only used for find positive or negitive chains."""
    chains_info = {}
    with open(gtf_file, "r") as f:
        for line in f:
            line = line.split("\t")
            chain = line[6]
            trans_name = line[8].split()[3].replace("\"","").split(".")[0]
            if trans_name not in chains_info:
                chains_info[trans_name] = chain
    return chains_info

def deal_indel(indel_file):
    adict = {}
    with open(indel_file, "r") as f:
        for line in f:
            line = line.rstrip("\n")
            # NM_004064  insertion c.356_357insTG  TG  2
            NM_id, flag, pos_info, nbase, length = line.split("\t")
            pos = pos_info.split("_")[0]
            if "insert" in flag:
                adict.setdefault(NM_id,[]).append(("insertion", pos, nbase ))
            elif "dup" in flag:
                adict.setdefault(NM_id,[]).append(("dup", pos, length))
            else:
                sys.exit("no type!")
    return adict

def dna_complement(sequence):
    ''' this is function docstring'''
    sequence = sequence.upper()
    sequence = sequence.replace('A', 't')
    sequence = sequence.replace('T', 'a')
    sequence = sequence.replace('C', 'g')
    sequence = sequence.replace('G', 'c')
    return sequence.upper()


def dna_reverse(sequence):
    ''' this is function docstring'''
    sequence = sequence.upper()
    return sequence[::-1]

def main():
    """docstring for main"""
    parser = _argparse()
    CDS_file = parser.cds_file
    gtf_file = parser.gtf_file
    output = open(parser.output_file, "w")

    cds_dict = deal_fasta(CDS_file)
    chains_info = deal_gtf(gtf_file)
    if parser.indel_file:
        insert_info = parser.indel_file
        indel_dict = deal_indel(insert_info)
        for each in cds_dict:
            chain = chains_info[each]
            cds_seq = cds_dict[each]
            if each in indel_dict:
                indel_info_list = indel_dict[each]
                sorted_list = sorted(indel_info_list, key=lambda item:item[1], reverse=True)
                for item in sorted_list:    
                    flag, pos, x = item
                    if chain == "-":
                        if flag == "insertion":
                            pos = int(re.match(r'c\.(\d+)\w+', pos).group(1)) # 'c.549dup' 
                            newseq = cds_seq[:pos-1]+ dna_reverse(dna_complement(x)) + cds_seq[pos-1:]
                        elif flag == "dup":
                            print("wow!!!!!!!!!!!!!!!!!!!!!!!!!!!!")
                            length = int(x)
                            pos = int(re.match(r'c\.(\d+)\-\w+', pos).group(1))  #  c.2042-747_*424dup
                            if cds_seq[pos-1].islower():
                                lowerseq = re.sub(r'[A-Z]','',cds_seq[pos-1:])  # replace upercase word.
                                newseq = cds_seq[:pos+length] + lowerseq + cds_seq[pos+length:]
                            else:
                                sys.exit("dup with wrong pos")
                        cds_dict[each] = newseq
                        newseq = cds_dict[each].replace("N", "")
                        output.write(">{}\n{}\n".format(each, newseq))
                    else:
                        if flag == "insertion":
                            pos = int(re.match(r'c\.(\d+)\w+', pos).group(1)) # 'c.549dup'
                            newseq = cds_seq[:pos]+ x + cds_seq[pos:]
                        elif flag == "dup":
                            length = int(x)
                            pos = int(re.match(r'c\.(\d+)\-\w+', pos).group(1))  #  c.2042-747_*424dup
                            if cds_seq[pos-1].islower():
                                lowerseq = re.sub(r'[A-Z]','',cds_seq[pos-1:])  # replace upercase word.
                                newseq = cds_seq[:pos+length] + lowerseq + cds_seq[pos+length:]
                            else:
                                sys.exit("dup with wrong pos")
                        cds_dict[each] = newseq
                        newseq = cds_dict[each].replace("N", "")
                        output.write(">{}\n{}\n".format(each, newseq))

            else:
                seq = cds_seq.replace("N", "")
                output.write(">{}\n{}\n".format(each,seq))
    else:
        for each in cds_dict:
            seq = cds_dict[each].replace("N", "")
            output.write(">{}\n{}\n".format(each, seq))
    output.close()


if __name__ == '__main__':
    main()
