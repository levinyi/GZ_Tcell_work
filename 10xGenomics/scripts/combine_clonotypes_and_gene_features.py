import sys,os


def usage():
    print("""
    usage:
        python {0} [file1] [file2] [file3]
    example:
        python {0} all_contig_annotations.csv clonotypes_add_consensus_annotations.csv Gene_expressed.count.csv >xxx.csv

    Update:
    20200529: created.
    """.format(os.path.basename(sys.argv[0])))

def deal_contig_file(contig_file):
    adict = {}
    with open(contig_file, "r") as f:
        for line in f:
            line = line.rstrip("\n")
            if line.startswith("barcode"):
                continue
            a = line.split(",")
            clonotype = a[-2]
            cell_barcode = a[0]
            if clonotype not in adict:
                adict.setdefault(clonotype,[]).append(cell_barcode)
            else:
                if cell_barcode not in adict.values():
                    adict.setdefault(clonotype, []).append(cell_barcode)
    return adict

def addtwodimdict(thedict, key_a, key_b, val):
    ''' this is a function to add two dimetion dict '''
    if key_a in thedict:
        thedict[key_a].update({key_b: val})
    else:
        thedict.update({key_a: {key_b: val}})
    return thedict

def deal_Gene_expressed(gene_expressed_file):
    express_dict = {}
    with open(gene_expressed_file, "r") as f:
        for line in f:
            line = line.rstrip("\n")
            if line.startswith(","):
                continue
            barcode,CD8A,CD8B,PDCD1,CD4 = line.split(",")
            CD8 = int(CD8A) + int(CD8B)
            addtwodimdict(express_dict, barcode, 'CD8',CD8)
            addtwodimdict(express_dict, barcode, 'PDCD1', int(PDCD1))
            addtwodimdict(express_dict, barcode, 'CD4', int(CD4))
    return express_dict


def main():
    if len(sys.argv) != 4:
        usage()
        sys.exit()

    contig_file = sys.argv[1]
    clonotype_file = sys.argv[2]
    gene_expressed_file = sys.argv[3]
    
    clone2barcode_dict = deal_contig_file(contig_file)
    express_dict = deal_Gene_expressed(gene_expressed_file)
    
    with open(clonotype_file, "r") as f:
        for line in f:
            line = line.rstrip("\n")
            if line.startswith("raw_clonotype_id"):
                continue
            a = line.split(",")
            clonotype = a[0]
            cell_barcodes = clone2barcode_dict[clonotype]
            CD8_max = 0
            CD4_max = 0
            PDCD1_max = 0
            for each in cell_barcodes:
                if express_dict[each]['CD8'] > CD8_max:
                    CD8_max = express_dict[each]['CD8']
                if express_dict[each]['CD4'] > CD4_max:
                    CD8_max = express_dict[each]['CD4']
                if express_dict[each]['PDCD1'] > PDCD1_max:
                    CD8_max = express_dict[each]['PDCD1']
            print(line,CD4_max,CD8_max,PDCD1_max)

if __name__ == '__main__':
    main()
