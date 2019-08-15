""" the """
import sys


def deal_blastn(blastn_result):
    """docstring for deal_blastn"""
    adict = {}
    with open(blastn_result, "r") as blastn_file:
        for line in blastn_file:
            if line.startswith("#"):
                continue
            line = line.rstrip("\n")
            record = line.split("\t")
            adict[record[0]] = record[1]
    return adict


def main():
    """docstring for main"""

    tra_file = sys.argv[1]
    trb_file = sys.argv[2]
    out_file = sys.argv[3]

    tra_dict = deal_blastn(tra_file)
    trb_dict = deal_blastn(trb_file)

    common_id = set(tra_dict.keys()) | set(trb_dict.keys())
    with open(out_file, "w") as out_put:
        for each in common_id:
            out_put.write(
                "{0}\t{1}\t{2}\n".format(
                    each,
                    tra_dict.get(each, 'NULL'),
                    trb_dict.get(each, "NULL")
                )
            )


if __name__ == '__main__':
    main()
