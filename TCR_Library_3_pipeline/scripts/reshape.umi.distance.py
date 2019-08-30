import sys
import re


def usage():
    """docstring for usage

    python $0 G30E1L2.merge.umi.with.distance.xls >G30E1L2.clonotype.umi.reads.xls

    2019-08-27: fix bug: enable find TRDV. c = re.compile(r'(TR[AB]V[0-9]+-?[0-9]?(DV\d)?)(.*)(TR[ABD]J.*)')
    2019-08-27: allowed TRDV/TRDJ in result.

    """


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


def main():
    """docstring for main"""
    TRV_gene = ["TRAV1-1", "TRAV1-2", "TRAV10", "TRAV11", "TRAV12-1", "TRAV12-2", "TRAV12-3", "TRAV13-1", "TRAV13-2", "TRAV14DV4", "TRAV16", "TRAV17", "TRAV18", "TRAV19", "TRAV2", "TRAV20", "TRAV21", "TRAV22", "TRAV23DV6", "TRAV24", "TRAV25", "TRAV26-1", "TRAV26-2", "TRAV27", "TRAV29DV5", "TRAV3", "TRAV30", "TRAV34", "TRAV35", "TRAV36DV7", "TRAV38-1", "TRAV38-2DV8", "TRAV39", "TRAV4", "TRAV40", "TRAV41", "TRAV5", "TRAV6", "TRAV7", "TRAV8-1", "TRAV8-2", "TRAV8-3", "TRAV8-4", "TRAV8-6", "TRAV8-7", "TRAV9-1", "TRAV9-2", "TRBV1", "TRBV10-1", "TRBV10-2", "TRBV10-3", "TRBV11-1", "TRBV11-2", "TRBV11-3", "TRBV12-1", "TRBV12-2", "TRBV12-3", "TRBV12-4",
                "TRBV12-5", "TRBV13", "TRBV14", "TRBV15", "TRBV16", "TRBV17", "TRBV18", "TRBV19", "TRBV2", "TRBV20-1", "TRBV21-1", "TRBV23-1", "TRBV24-1", "TRBV25-1", "TRBV26", "TRBV27", "TRBV28", "TRBV29-1", "TRBV3-1", "TRBV3-2", "TRBV30", "TRBV4-1", "TRBV4-2", "TRBV4-3", "TRBV5-1", "TRBV5-3", "TRBV5-4", "TRBV5-5", "TRBV5-6", "TRBV5-7", "TRBV5-8", "TRBV6-1", "TRBV6-2", "TRBV6-3", "TRBV6-4", "TRBV6-5", "TRBV6-6", "TRBV6-7", "TRBV6-8", "TRBV6-9", "TRBV7-1", "TRBV7-2", "TRBV7-3", "TRBV7-4", "TRBV7-6", "TRBV7-7", "TRBV7-8", "TRBV7-9", "TRBV9", "TRAV8-5", "TRAV28", "TRBV7-5", "TRBV5-2", "TRAV33", "TRBV22-1"]  # last six genes are pseudo genes. added by hand.
    TRJ_gene = ["TRAJ1", "TRAJ10", "TRAJ11", "TRAJ12", "TRAJ13", "TRAJ14", "TRAJ15", "TRAJ16", "TRAJ17", "TRAJ18", "TRAJ19", "TRAJ2", "TRAJ20", "TRAJ21", "TRAJ22", "TRAJ23", "TRAJ24", "TRAJ25", "TRAJ26", "TRAJ27", "TRAJ28", "TRAJ29", "TRAJ3", "TRAJ30", "TRAJ31", "TRAJ32", "TRAJ33", "TRAJ34", "TRAJ35", "TRAJ36", "TRAJ37", "TRAJ38", "TRAJ39", "TRAJ4", "TRAJ40", "TRAJ41", "TRAJ42", "TRAJ43", "TRAJ44", "TRAJ45", "TRAJ46", "TRAJ47", "TRAJ48", "TRAJ49", "TRAJ5", "TRAJ50", "TRAJ51", "TRAJ52", "TRAJ53", "TRAJ54", "TRAJ55", "TRAJ56", "TRAJ57", "TRAJ58", "TRAJ59", "TRAJ6", "TRAJ60", "TRAJ61", "TRAJ7", "TRAJ8", "TRAJ9", "TRBJ1-1", "TRBJ1-2", "TRBJ1-3", "TRBJ1-4", "TRBJ1-5", "TRBJ1-6", "TRBJ2-1", "TRBJ2-2", "TRBJ2-2P", "TRBJ2-3", "TRBJ2-4", "TRBJ2-5", "TRBJ2-6", "TRBJ2-7"]
    c = re.compile(r'(TR[ABD]V[0-9]+-?[0-9]?(DV\d)?)(.*)(TR[ABD]J.*)')  # fix bugs *TR[AB]V*

    adict = {}
    with open(sys.argv[1], "r") as f:
        for line in f:
            line = line.rstrip("\n")
            if line.startswith("uniqueCloneType"):
                continue
            clonotype, umi, readsNumber = line.split("\t")
            two_dim_dict(adict, clonotype, 'uminumber', 1)
            two_dim_dict(adict, clonotype, 'readsNumber', int(readsNumber))
    print("Clonotype\tTRV\tCDR3\tTRJ\tUMIcount\tReadsNumber")
    for clonotype in adict:
        p = c.match(clonotype)
        if p:
            v_gene = p.group(1)
            cdr3 = p.group(3)
            J_gene = p.group(4)

            if v_gene in TRV_gene and J_gene in TRJ_gene:
                print("%s\t%s\t%s\t%s\t%s\t%s" % (clonotype, v_gene, cdr3, J_gene, adict[clonotype]['uminumber'], adict[clonotype]['readsNumber']))
            elif v_gene in TRV_gene and J_gene not in TRJ_gene:
                if 'TRD' in v_gene or 'TRD' in J_gene:
                    print("%s\t%s\t%s\t%s\t%s\t%s" % (clonotype, v_gene, cdr3, J_gene, adict[clonotype]['uminumber'], adict[clonotype]['readsNumber']))
                else:
                    print("%s\t%s\t%s\t%s\t%s\t%s" % (clonotype, v_gene, cdr3, J_gene, adict[clonotype]['uminumber'], adict[clonotype]['readsNumber']))
                    # sys.exit("NULL_J")
                    sys.stderr.write("This TRJ not in your dict: {}\n".format(clonotype))
            elif v_gene not in TRV_gene and J_gene in TRJ_gene:
                if 'TRD' in v_gene or 'TRD' in J_gene:
                    print("%s\t%s\t%s\t%s\t%s\t%s" % (clonotype, v_gene, cdr3, J_gene, adict[clonotype]['uminumber'], adict[clonotype]['readsNumber']))
                else:
                    print("%s\t%s\t%s\t%s\t%s\t%s" % (clonotype, v_gene, cdr3, J_gene, adict[clonotype]['uminumber'], adict[clonotype]['readsNumber']))
                    # sys.exit("NULL_V")
                    sys.stderr.write("This TRV not in your dict: {}\n".format(clonotype))
            elif v_gene not in TRV_gene and J_gene not in TRJ_gene:
                if 'TRD' in v_gene or 'TRD' in J_gene:
                    print("%s\t%s\t%s\t%s\t%s\t%s" % (clonotype, v_gene, cdr3, J_gene, adict[clonotype]['uminumber'], adict[clonotype]['readsNumber']))
                else:
                    print("%s\t%s\t%s\t%s\t%s\t%s" % (clonotype, v_gene, cdr3, J_gene, adict[clonotype]['uminumber'], adict[clonotype]['readsNumber']))
                    sys.stderr.write("This TRV and TRJ are not in your dict: {}\n".format(clonotype))


if __name__ == '__main__':
    main()
