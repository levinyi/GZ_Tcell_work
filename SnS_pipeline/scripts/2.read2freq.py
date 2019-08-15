''' this is the module docstring'''
import sys


def two_dim_dict(thedict, key_a, key_b, value):
    ''' this is '''
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

    input_file = sys.argv[1]
    output_file = sys.argv[2]

    a_dict = {}
    with open(input_file, "r") as input_file:
        for line in input_file:
            if 'NULL' in line:
                continue
            line = line.rstrip("\n")
            record = line.split("\t")
            two_dim_dict(a_dict, record[1], record[2], 1)

    total_pairs = 0
    perfect_match = 0
    mispairing_match = 0
    with open(output_file, "w") as op:
        for alpha in a_dict:
            for beta in a_dict[alpha]:
                if alpha == beta:
                    perfect_match += 1
                else:
                    mispairing_match += 1
                total_pairs += 1
                op.write("{}\t{}\t{}\n".format(alpha, beta, a_dict[alpha][beta]))
    print("{0}\t{1}\t{2}".format(total_pairs, perfect_match, mispairing_match))

if __name__ == '__main__':
    main()
