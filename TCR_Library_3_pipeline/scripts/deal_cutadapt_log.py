import sys
import os
import argparse


def _argparse():
    parser = argparse.ArgumentParser(description="This is description")
    parser.add_argument('-l', '--log', action='store', dest='log_file', required=True, help="input read1 file")
    parser.add_argument('-d', '--dir', action='store', dest='result_dir', default='./', help="input read2 file")
    parser.add_argument('-f',  '--flag', action='store_true', dest='flag', default=False, help="means to contains -l flag in output.")

    parser.add_argument('-v',  '--version', action='version', version='%(prog)s 0.1')
    return parser.parse_args()


def print_current_time():
    time_stamp = datetime.datetime.now()
    return time_stamp.strftime('%Y.%m.%d-%H:%M:%S')


def main():
    parser = _argparse()

    logfile = parser.log_file
    result_dir = parser.result_dir
    # print parser.flag
    if parser.flag == True:
        print parser.flag
    else:
        print parser.flag
    file_name = os.path.basename(logfile).split(".")[0]

    length_list = []
    with open(logfile, "r") as f, open(result_dir + '/' + file_name + ".p1.fq", "w") as OUT1, open(result_dir + '/' + file_name + '.p2.fq', "w") as OUT2:
        for line in f:
            line = line.rstrip("\n")
            c = line.split("\t")
            length_list.append(len(c))
            if len(c) == 4:  # no adapter in line
                if parser.flag == True: # if flag=true,print to output.
                    OUT1.write("@%s\n%s\n+\n%s\n" % (c[0], c[2], c[3]))
            elif len(c) == 11:  # with adapter line
                if len(c[4]) == 0:  # adapter starts from front.
                    OUT2.write("@%s\n%s\n+\n%s\n" % (c[0], c[6], c[10]))
                elif len(c[6]) == 0:  # adapter is in the end.
                    OUT1.write("@%s\n%s\n+\n%s\n" % (c[0], c[4], c[8]))
                else:  # adapter is in the middle.
                    OUT1.write("@%s\n%s\n+\n%s\n" % (c[0], c[4], c[8]))
                    OUT2.write("@%s\n%s\n+\n%s\n" % (c[0], c[6], c[10]))
            else:
                print("Attention: your log file contains : %s columns." % len(c))
    print("finished")
    
if __name__ == '__main__':
    main()
