import sys


contig_file = sys.argv[1]
       

class ClassName(object):
    """docstring for ClassName"""
    def __init__(self, cell_barcode,vdj_A,vdj_B,umi_A,umi_B):
        self.cell_barcode = cell_barcode
        self.vdj_A = vdj_A
        self.vdj_B = vdj_B
        self.umi_A = umi_A
        self.umi_B = umi_B



def deal_file():
    with open() as f:
        for line in f:
            line = line.