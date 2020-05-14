import sys


input_file = sys.argv[1]
# OM_ID = sys.argv[1]

# print("OM_ID\tGene\tTranscript_Ref\tClassification\tAlternation_type1\tAlternationtype2\tChromosome\tPositionStart\tPositionEnd\tGenomicDNAchange\tCodingDNAchange\tAAchange\tVaf\tMeanCoverage")
with open(input_file, "r") as f:
    for line in f:
        if line.startswith("CHROM"):
            continue
        chrom = line.split("\t")[0]
        pos = line.split("\t")[1]
        ExonicFunc = line.split("\t")[11]
        AAChange = line.split("\t")[12]
        REF = line.split("\t")[3]
        ALT = line.split("\t")[4]
        if ExonicFunc != "." and AAChange != "." and ExonicFunc != "synonymous SNV" and AAChange != "UNKNOWN":
            print("chr{}\t{}\t{}\t{}\t{}\t{}\t{}".format(chrom, pos, pos, REF, ALT, ExonicFunc, "\t".join(AAChange.split(",")[0].split(":"))))



