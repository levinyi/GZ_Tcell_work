# hla_scan -b NA12155.chr6.bam -d HLA-ALL.IMGT -v 38
hla_scan -t 60 -l CP60005015-C01-WES_R1.fastq -r CP60005015-C01-WES_R2.fastq -d /cygene/database/imgt-hla/db/HLA-ALL.IMGT -g HLA-A HLA-B HLA-C HLA-E HLA-F HLA-G MICA MICB HLA-DMA HLA-DMB HLA-DOA HLA-DOB HLA-DPA1 HLA-DPB1 HLA-DQA1 HLA-DQB1 HLA-DRA HLA-DRB1 HLA-DRB5 TAP1 TAP2

#time-consummer

for i in {"HLA-A","HLA-B","HLA-C","HLA-E","HLA-F","HLA-G","MICA","MICB","HLA-DMA","HLA-DMB","HLA-DOA","HLA-DOB","HLA-DPA1","HLA-DPB1","HLA-DQA1","HLA-DQB1","HLA-DRA","HLA-DRB1","HLA-DRB5","TAP1","TAP2"};do echo "nohup hla_scan -t 10 -l CP60005015-C01-WES_R1.fastq -r CP60005015-C01-WES_R2.fastq -d /cygene/database/imgt-hla/db/HLA-ALL.IMGT -g $i >$i.out.txt &" ;done

nohup hla_scan -t 60 -l CP60005015-C01-WES_R1.fastq -r CP60005015-C01-WES_R2.fastq -d /cygene/database/imgt-hla/db/HLA-ALL.IMGT -g HLA-A >HLA-A.out.txt &
nohup hla_scan -t 60 -l CP60005015-C01-WES_R1.fastq -r CP60005015-C01-WES_R2.fastq -d /cygene/database/imgt-hla/db/HLA-ALL.IMGT -g HLA-B >HLA-B.out.txt &
nohup hla_scan -t 60 -l CP60005015-C01-WES_R1.fastq -r CP60005015-C01-WES_R2.fastq -d /cygene/database/imgt-hla/db/HLA-ALL.IMGT -g HLA-C >HLA-C.out.txt &
nohup hla_scan -t 60 -l CP60005015-C01-WES_R1.fastq -r CP60005015-C01-WES_R2.fastq -d /cygene/database/imgt-hla/db/HLA-ALL.IMGT -g HLA-E >HLA-E.out.txt &
nohup hla_scan -t 60 -l CP60005015-C01-WES_R1.fastq -r CP60005015-C01-WES_R2.fastq -d /cygene/database/imgt-hla/db/HLA-ALL.IMGT -g HLA-F >HLA-F.out.txt &
nohup hla_scan -t 60 -l CP60005015-C01-WES_R1.fastq -r CP60005015-C01-WES_R2.fastq -d /cygene/database/imgt-hla/db/HLA-ALL.IMGT -g HLA-G >HLA-G.out.txt &
nohup hla_scan -t 60 -l CP60005015-C01-WES_R1.fastq -r CP60005015-C01-WES_R2.fastq -d /cygene/database/imgt-hla/db/HLA-ALL.IMGT -g MICA >MICA.out.txt &
nohup hla_scan -t 60 -l CP60005015-C01-WES_R1.fastq -r CP60005015-C01-WES_R2.fastq -d /cygene/database/imgt-hla/db/HLA-ALL.IMGT -g MICB >MICB.out.txt &
nohup hla_scan -t 60 -l CP60005015-C01-WES_R1.fastq -r CP60005015-C01-WES_R2.fastq -d /cygene/database/imgt-hla/db/HLA-ALL.IMGT -g HLA-DMA >HLA-DMA.out.txt &
nohup hla_scan -t 60 -l CP60005015-C01-WES_R1.fastq -r CP60005015-C01-WES_R2.fastq -d /cygene/database/imgt-hla/db/HLA-ALL.IMGT -g HLA-DMB >HLA-DMB.out.txt &
nohup hla_scan -t 60 -l CP60005015-C01-WES_R1.fastq -r CP60005015-C01-WES_R2.fastq -d /cygene/database/imgt-hla/db/HLA-ALL.IMGT -g HLA-DOA >HLA-DOA.out.txt &
nohup hla_scan -t 60 -l CP60005015-C01-WES_R1.fastq -r CP60005015-C01-WES_R2.fastq -d /cygene/database/imgt-hla/db/HLA-ALL.IMGT -g HLA-DOB >HLA-DOB.out.txt &
nohup hla_scan -t 60 -l CP60005015-C01-WES_R1.fastq -r CP60005015-C01-WES_R2.fastq -d /cygene/database/imgt-hla/db/HLA-ALL.IMGT -g HLA-DPA1 >HLA-DPA1.out.txt &
nohup hla_scan -t 60 -l CP60005015-C01-WES_R1.fastq -r CP60005015-C01-WES_R2.fastq -d /cygene/database/imgt-hla/db/HLA-ALL.IMGT -g HLA-DPB1 >HLA-DPB1.out.txt &
nohup hla_scan -t 60 -l CP60005015-C01-WES_R1.fastq -r CP60005015-C01-WES_R2.fastq -d /cygene/database/imgt-hla/db/HLA-ALL.IMGT -g HLA-DQA1 >HLA-DQA1.out.txt &
nohup hla_scan -t 60 -l CP60005015-C01-WES_R1.fastq -r CP60005015-C01-WES_R2.fastq -d /cygene/database/imgt-hla/db/HLA-ALL.IMGT -g HLA-DQB1 >HLA-DQB1.out.txt &
nohup hla_scan -t 60 -l CP60005015-C01-WES_R1.fastq -r CP60005015-C01-WES_R2.fastq -d /cygene/database/imgt-hla/db/HLA-ALL.IMGT -g HLA-DRA >HLA-DRA.out.txt &
nohup hla_scan -t 60 -l CP60005015-C01-WES_R1.fastq -r CP60005015-C01-WES_R2.fastq -d /cygene/database/imgt-hla/db/HLA-ALL.IMGT -g HLA-DRB1 >HLA-DRB1.out.txt &
nohup hla_scan -t 60 -l CP60005015-C01-WES_R1.fastq -r CP60005015-C01-WES_R2.fastq -d /cygene/database/imgt-hla/db/HLA-ALL.IMGT -g HLA-DRB5 >HLA-DRB5.out.txt &
nohup hla_scan -t 60 -l CP60005015-C01-WES_R1.fastq -r CP60005015-C01-WES_R2.fastq -d /cygene/database/imgt-hla/db/HLA-ALL.IMGT -g TAP1 >TAP1.out.txt &
nohup hla_scan -t 60 -l CP60005015-C01-WES_R1.fastq -r CP60005015-C01-WES_R2.fastq -d /cygene/database/imgt-hla/db/HLA-ALL.IMGT -g TAP2 >TAP2.out.txt &

