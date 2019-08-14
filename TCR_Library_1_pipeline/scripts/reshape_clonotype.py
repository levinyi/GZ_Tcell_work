import os 
import re

def split_clonotype(clonotype):
    c = re.compile(r'(TR[AB]V[0-9]+-?[0-9]?(DV\d)?)(.*)(TR[ABD]J.*)')
    p = c.match(clonotype)
    if p:
        v_gene = p.group(1)
        cdr3 = p.group(3)
        J_gene = p.group(4)
        return v_gene,cdr3,J_gene

