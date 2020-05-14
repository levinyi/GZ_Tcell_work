import sys
from Bio import SeqIO


def deal_genome(genome_file):
    genome_dict = {}
    for record in SeqIO.parse(genome_file, "fasta"):
        genome_dict[str(record.id)] = str(record.seq).upper()
    return genome_dict


def main():
    genome_file = sys.argv[1]
    STF_file = sys.argv[2]
    indel_output = sys.argv[3]
    genome_dict = deal_genome(genome_file)
    with open(STF_file, "r") as f, open(indel_output, "w") as f2:
        for line in f:
            line = line.rstrip("\n")
            if line.startswith("OM_ID"):
                continue
            chrom, pos1, pos2, ref, alt, anno, gene, transref, exon, Coding_DNA_change, aachange = line.split("\t")
            pos = int(pos1)
            # chrom = 'chr'+ chrom
            old_seq = genome_dict[chrom][ pos - 10 : pos + 10 ]
            
            if genome_dict[chrom][pos-1] != ref[0]:
                print("REF is not Correct,please check! {}\t{}".format(ref[0],genome_dict[chrom][pos-1]),end="\t")
            
            if anno.endswith(("SNV","stoploss")):
                new_seq = genome_dict[chrom][:pos-1] + alt +genome_dict[chrom][pos:]
            elif anno.endswith("stopgain"):
                if len(ref) == len(alt): # snp
                    new_seq = genome_dict[chrom][:pos-1] + alt +genome_dict[chrom][pos:]
                elif len(ref) <len(alt): # insertion
                    new_seq = genome_dict[chrom][:pos-1] + genome_dict[chrom][pos-1].lower() + genome_dict[chrom][pos:]
                    f2.write("{}\tinsertion\t{}\t{}\t{}\n".format(transref,Coding_DNA_change,alt[1:],len(alt[1:])))
            elif anno.endswith("insertion"):
                new_seq = genome_dict[chrom][:pos-1] + genome_dict[chrom][pos-1].lower() + genome_dict[chrom][pos:]
                f2.write("{}\tinsertion\t{}\t{}\t{}\n".format(transref,Coding_DNA_change,alt[1:],len(alt[1:])))
            
            elif anno.endswith("deletion"):
                new_seq = genome_dict[chrom][:pos] + (len(ref)-1)*"N" + genome_dict[chrom][pos+len(ref)-1:]
            
            else:
                sys.exit("no such annotation info had seen before.!!!!!")
            testseq = new_seq[pos-10: pos+10]
            print(line, old_seq, testseq)
            genome_dict[chrom] = new_seq
    
    with open("genome_replaced_mutation.fa", "w") as output:
        for each in genome_dict:
            output.write(">{}\n{}\n".format(each, genome_dict[each]))


if __name__ == '__main__':
    main()