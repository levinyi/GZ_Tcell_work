import sys
import os
import json


def addtwodimdict(thedict, key_a, key_b, val):
    ''' this is a function to add two dimetion dict '''
    if key_a in thedict:
        thedict[key_a].update({key_b: val})
    else:
        thedict.update({key_a: {key_b: val}})
    return thedict


def deal_STAR_log_file():
    mapping_dict = {}
    for root, dirs, files in os.walk("."):
        for each_file in files:
            if each_file.endswith("Log.final.out"):
                # print(each_file)
                name = each_file.split(".")[0]
                sub_name = ".".join(each_file.split(".")[1:-3])
                # print(name, sub_name)
                sample_dict = {}
                with open(each_file, "r") as f:
                    for line in f:
                        line = line.rstrip("\n")
                        if "Number of input reads" in line:
                            total_reads = line.split("|")[1].lstrip("\t")
                        if "Uniquely mapped reads number" in line:
                            mapped_reads = line.split("|")[1].lstrip("\t")
                        if "Uniquely mapped reads %" in line:
                            mapped_rate = line.split("|")[1].lstrip("\t")
                addtwodimdict(sample_dict, name, "total_reads", total_reads)
                addtwodimdict(sample_dict, name, "mapped_reads", mapped_reads)
                addtwodimdict(sample_dict, name, "mapped_rate", mapped_rate)
                addtwodimdict(mapping_dict, name, sub_name, sample_dict[name])
    #print(json.dumps(mapping_dict, indent=4))
    return mapping_dict


def deal_ribosomal():
    '''scan work directory to find *bbduk.out file as input.
        read 'Input: 12896786 reads'
        'Contaminants: 360488 reads (2.80%)' as output.
    '''
    ribosomal_dict = {}
    for root, dirs, files in os.walk("."):
        for each_file in files:
            if each_file.endswith(".bbduk.out"):
                file_name = ".".join(each_file.split(".")[:-2])
                #print("find file: {}\t basename: {}".format(each_file, file_name))
                with open(each_file,"r") as f:
                    for line in f:
                        line = line.rstrip("\n")
                        if line.startswith("Input"):
                            # this is the Input : xxx reads line.
                            input_reads = line.split()[1]
                        if line.startswith("Contaminants"):
                            # find Contaminants line
                            contamin_reads = line.split()[1]
                            contamin_rate = line.split()[3].strip("(").strip(")")
                    addtwodimdict(ribosomal_dict, file_name, 'input_reads', input_reads)
                    addtwodimdict(ribosomal_dict, file_name, 'contamin_reads', contamin_reads)
                    addtwodimdict(ribosomal_dict, file_name, 'contamin_rate', contamin_rate)
    #print(json.dumps(ribosomal_dict, indent=4))
    return(ribosomal_dict)


def main():
    mapped_dict = deal_STAR_log_file()
    ribosomal_dict = deal_ribosomal()
    # write to file
    # header 
    print("Sample\tTotal_Reads\tmapped_mm10_reads\trate\tmm10.unmapped.to.hg38.total.reads\tmapped_reads\trate\ttotal.reads.for.hg38\tmapped_reads\trate\ttotal.reads.for.EBV\tmapped_reads\trate\ttotal_reads.for.ribosomal\tcontamin_reads\tcontamin_rate")
    header = ["mm10","mm10.unmapped.to.hg38","hg38", "EBV"]
    for name, v in mapped_dict.items():
        output_content = [name,]
        for each in header:
            output_content.append(v[each].get("total_reads"))
            output_content.append(v[each].get("mapped_reads"))
            output_content.append(v[each].get("mapped_rate"))
        output_content.append(ribosomal_dict[name]['input_reads'])
        output_content.append(ribosomal_dict[name]['contamin_reads'])
        output_content.append(ribosomal_dict[name]['contamin_rate'])
        print("\t".join(output_content))

if __name__ == '__main__':
    main()

