import sys
sys.path.append("/cygene/script/python_packages")
import edit_distance


def usage():
    """docstring for usage
    python .py /cygene/work/00.test/pipeline/TCR_Library_1_pipeline/database/VD
    Jdb-WT1.xls G36E1L2.clonotype.umi.reads.xls 
    """

db_file = sys.argv[1]
clonotype_file = sys.argv[2]


class DB_WT1():
    '''def __init__(self,Gene, CDR3, V, J, Species,MHC_A,MHC_B,MHC_class,Epitope,Epitope_gene,Epitope_species,Score): '''
    def __init__(self, Gene, CDR3, V, J, Species, Score): 
        """docstring for __init__"""
        self.Gene = Gene
        self.CDR3 = CDR3
        self.V = V
        self.J = J
        self.Species = Species
        self.Score = Score

    def show_detail(self):
        return self.Gene, self.CDR3, self.V, self.J

CDR3_db_list = []
with open(db_file, "r") as f:
    for line in f:
        line = line.rstrip("\n")
        if line.startswith("Gene"):
            continue
        Gene, CDR3, V, J, Species, MHC_A, MHC_B, MHC_class, Epitope, Epitope_gene, Epitope_species, Score = line.split("\t")
        aninstance = DB_WT1(Gene, CDR3, V, J, Species, Score)
        CDR3_db_list.append(aninstance)


with open(clonotype_file, "r") as f2:
    for line in f2:
        line = line.rstrip("\n")
        if line.startswith("Clonotype"):
            continue
        Clonotype, TRV, CDR3, TRJ, UMIcount, ReadsNumber = line.split("\t")
        closed_cdr3 = {}
        for each in CDR3_db_list:
            l = edit_distance.minEditDist(CDR3, each.CDR3)
            if l <= 2:
                closed_cdr3[each] = l

        if len(closed_cdr3) > 0:
            print line, [each.show_detail() for each in closed_cdr3]
