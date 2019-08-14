import sys
import os


def usage():
    '''
    usage: python M1_relaxed.py <M1 file> > M1_relaxed.xls
    example: python M1_relaxed.py /cygene/work/27.G13/G13AB/result/G13E3L8.pairs.freq.M1.xls > /cygene/work/27.G13/G13AB/result/G13E3L8.pairs.freq.relaxed.M1.xls
    
    2019-06-20: fix some bugs.
    '''
    pass


def compare(c, d):
    # c = a.split("\t")
    # d = b.split("\t")
    if c[0] == d[0]:
        Alpha = c[0]
        if int(c[10]) >= int(d[10]):
            # Beta = c[10] bug
            Beta = c[4]
        else:
            # Beta = d[10] bug
            Beta = d[4]
        ReadCount = int(c[8]) + int(d[8])
        A_total_reads = int(c[9])
        B_total_reads = int(c[10]) + int(d[10])
    elif c[4] == d[4]:
        if int(c[9]) >= int(d[9]):
            Alpha = c[0]
        else:
            Alpha = d[0]
        Beta = c[4]
        ReadCount = int(c[8]) + int(d[8])
        A_total_reads = int(c[9]) + int(d[9])
        B_total_reads = int(c[10])
    else:
        if int(c[9]) >= int(d[9]):
            Alpha = c[0]
        else:
            Alpha = d[0]

        # if c[10] >= d[10]: bug
        if int(c[10]) >= int(d[10]):
            Beta = c[4]
        else:
            Beta = d[4]
        ReadCount = int(c[8]) + int(d[8])
        A_total_reads = int(c[9]) + int(d[9])
        B_total_reads = int(c[10]) + int(d[10])
        # print("attention")

    M1AScore = float(ReadCount) / A_total_reads
    M1BScore = float(ReadCount) / B_total_reads
    M1Score = float(M1AScore * M1BScore)
    merge_list = []
    merge_list.extend([Alpha, c[1], c[2], c[3], Beta, c[5], c[6], c[7], str(ReadCount), str(A_total_reads), str(B_total_reads), str(M1AScore), str(M1BScore), str(M1Score)])
    
    return merge_list


def main():
    adict = {}
    with open(sys.argv[1], "r") as f:
        for line in f:
            line = line.rstrip("\n")
            if line.startswith("Alpha"):
                continue
            Alpha, TRAV, CDR3Alpha, TRAJ, Beta, TRBV, CDR3Beta, TRBJ, ReadCount, A_total_reads, B_total_reads, M1AScore, M1BScore, M1Score = line.split("\t")
            l = line.split("\t")
            pairs = TRAV + CDR3Alpha + TRAJ + '_' + TRBV + CDR3Beta + TRBJ
            if pairs not in adict:
                adict[pairs] = l
            else:
                new_pairs = compare(adict[pairs], l)
                adict[pairs] = new_pairs
    for k in adict:
        print "\t".join(adict[k])

if __name__ == '__main__':
    main()
