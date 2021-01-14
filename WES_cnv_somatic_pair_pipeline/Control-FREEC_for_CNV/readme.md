download 
tar -zxvf FREEC-11.6.tar.gz
cd FREEC-11.6/src
make
sudo ln -s /cygene/software/biosoftware/FREEC-11.6/src/freec /usr/local/bin/

test:
mkdir Control-FREEC_for_CNV

cd Control-FREEC_for_CNV

## prepare data
cp /cygene/software/biosoftware/FREEC-11.6/data/config_exome.txt ./
ln -s xxx.bam ./mySample.bam
ln -s xxx.bam ./myControl.bam

####
# freec -conf <config file> -sample <mySample.bam> -control <myControl.bam>
freec -conf config_exome.txt
