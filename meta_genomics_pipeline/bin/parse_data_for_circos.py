import glob
import os
import pandas as pd

def make_gene_list():
    """
    parse out the gene names of interest
    """
    os.system('mkdir -p org_gff')
    ########### change final.contigs.long.fa to my file.
    os.system("grep '>' final.contigs.long.fa |cut -f1 -d ' '|cut -f2 -d '>' > org_list")
    ########### change metagG.gff to my file.
    os.system("for org in `cat org_list`; do grep -w $org metagG.gff > $org.subset.gff; done")
    ###########
    os.system('mv *subset.gff org_gff')

def loop_over_gff(directory):
    """
    loop over gff files within a directory
    """
    files=glob.glob(directory+'*.gff')
    gffHash={}
    for f in files:
        tmp=process_gff_file(f)
        gffHash[tmp['genome']]=tmp
    return(gffHash)


def process_gff_file(in_file):
    """
    process gff file line by line
    """
    handle=open(in_file)
    outDict={}
    outDict['genes']={}
    for line in handle:
        if line.startswith('#'):
            items=line.split(' ')
            outDict['genome']=items[1].strip()
            outDict['genome_start']=0
            outDict['genome_end']=items[-1].strip()
        elif line.startswith('>'):
            break
        else:
            g={}
            items=line.split('\t')
            g['start']=items[3]
            g['end']=items[4]
            g['orientation']=items[6]
            g['id']=items[8].split(';')[0].split('=')[1]
            g['def']=items[8].split(';')[-1].split('=')[1].strip()
            outDict['genes'][ g['id']]=g
    return outDict

def make_karyotype(gffHash, outfile):
    """
    make circos karyotype file
    """
    out=open(outfile, 'w')
    for i, genome in enumerate(gffHash):
        outlist=['chr', '-', genome, str(i+1),str(gffHash[genome]['genome_start']),
                 str(gffHash[genome]['genome_end']), 'chr'+str(i+1)]
        out.write(' '.join(outlist))
        out.write('\n')
    out.close()
    
    
    
def make_gene_tile(gffHash, outfile):    
    """
    make circos gene layout file colored based on orientation
    """
    outf=open(outfile+'forward.txt', 'w')
    outr=open(outfile+'reverse.txt', 'w')

    for i, genome in enumerate(gffHash):
        for gene in gffHash[genome]['genes']:
            g=gffHash[genome]['genes'][gene]
            o=g['orientation']
            if o=='-':
                ort=-1
                outlist=[genome, g['start'], g['end'], str(ort)]
                outr.write(' '.join(outlist))
                outr.write('\n')

            else:
                ort=1
                outlist=[genome, g['start'], g['end'], str(ort)]
                outf.write(' '.join(outlist))
                outf.write('\n')
    outf.close()
    outr.close()

def make_count_heatmap(count_file, gffHash):
    """
    given count files with gene name and count-- create hash for heatmap
    """
    df=pd.read_table(count_file, index_col='transcript')
    out=open(count_file+'.circos','w')
    for i, genome in enumerate(gffHash):
        for gene in gffHash[genome]['genes']:
            g=gffHash[genome]['genes'][gene]
            tpm=df.loc[gene,'count']
            outlist=[genome, g['start'], g['end'], str(tpm)]
            out.write(' '.join(outlist))
            out.write('\n')
    out.close()
                
                
def main():
    # modified by dushiyi
    name = sys.argv[1]
    make_gene_list()   
    gffHash=loop_over_gff('org_gff/')
    make_karyotype(gffHash, 'metag.karyotype.txt')
    make_gene_tile(gffHash, 'metag.genes.orientation.')
    for count in glob.glob('*counts'):
        make_count_heatmap(count, gffHash)
    print('done')

if __name__ == '__main__':
    main()
