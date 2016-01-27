from numpy import *
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord 
from Bio.Seq import Seq

# takes a pedline and processes it wrt a reference genome
# removes the first 6 columns
# collapses hetero snps and makes 2 versions of half size
# supressed has all hetero snps reverted to reference
# amplified has all hetero snps biased to variant
# also counts identities, homo, hetero snps and no-recalls 

def pedlineProcess(genome, index, pedline) :
    ncols = len(pedline)
    amplified = empty((ncols-6)/2, dtype=str )
    supressed = empty((ncols-6)/2, dtype=str )
    name = pedline[0]
    homo=0
    hetero=0
    identity=0
    norecall=0

    for i in range(6,ncols,2) :
        # number of snp after collapsing alelles
        nsnp=(i-6)/2

        # reference for this snp
        reference = genome[index[nsnp][0]] [index[nsnp][1]-1] 
        #print reference, pedline[i], pedline[i+1]
        #print name ,index[nsnp][0], index[nsnp][1], reference, pedline[i], pedline[i+1]

        # homozygous case
        if pedline[i] == pedline[i+1] :
            homo += 1
            amplified[nsnp] = pedline[i]
            supressed[nsnp] = pedline[i]

            if pedline[i] == reference:
                identity += 1
            elif pedline[i] == '0' :
                norecall += 1
              
        #heterozygous case
        else : 
            hetero += 1

            if pedline[i] == reference :
                amplified[nsnp] = pedline[i+1]
            elif pedline[i+1] == reference : 
                amplified[nsnp] = pedline[i]
            else :
                print 'Hetero but no allele is identical to reference. Wow, thats rare!'
                amplified[nsnp] = '0'

            supressed[nsnp] = reference;
            
    return name, supressed, amplified, (identity,homo,hetero,norecall)



# takes a (processed) line from a ped file 
# and mutates a genome accordingly 

def mutate (genome, index, pedline) :
    nsnps = len(pedline)
    newgenome=[0]

    ## ADD EXCEPTION/ERROR HERE!
     
    # make mutable copy of the genome 
    for i in range(1,13) :
        newgenome.append( genome[i].tomutable())

    # mutate the genome for each snp
    for i in range(nsnps) :
        chromosome = index[i][0]
        position   = index[i][1]-1
        mutateto   = pedline[i]
        
        if mutateto!='0' :
            newgenome[chromosome][position] = mutateto 
    for i in range(1,13) :
        newgenome[i] = SeqRecord (newgenome[i].toseq(), id='chr'+str(i), description='')

    return newgenome


def chrtoint(string) : return int(string[3:])

def parseGff(gffFilename) :
    gffint = loadtxt(gffFilename, dtype=int, usecols=(0,3,4), converters={0:chrtoint}) 
    # strand
    gffstr = loadtxt(gffFilename, dtype=str, usecols=([6]) )
    # attributes
    gffatt = loadtxt(gffFilename, dtype=str, usecols=([8]), comments=';' )
    tmp=[]
    for line in gffatt :
        tmp.append(line.split(':'))
    gffatt=array(tmp)
    tmp=0

    annotation={}

    for i in range(len(gffstr)) :
        if not annotation.has_key(gffatt[i][0]) :
            annotation[gffatt[i][0]] = [gffint[i][0], gffstr[i], (gffint[i][1], gffint[i][2])]
        else :
            annotation[gffatt[i][0]].append( (gffint[i][1],gffint[i][2]))

    return annotation

from Bio.Alphabet import generic_dna, generic_protein

#splits the genome to individual genes
def splitGenome(genome, annotation) :
    allgenes=[]
    for gene in annotation.keys()   :
        chromosome = annotation[gene][0]
        strand = annotation[gene][1]
        sequence=Seq('')
        for cds in annotation[gene][2:] :
            if strand == '+' :
                sequence += genome[chromosome].seq[cds[0]-1:cds[1]]
            elif strand == '-' :
                sequence +=  genome[chromosome].seq[cds[0]-1:cds[1]].reverse_complement()

        if sequence[:3] != 'ATG' :
            print gene, 'cds not starting with ATG!' , sequence

        allgenes.append( SeqRecord ( sequence.translate(),id=gene[3:]) )

    return allgenes  


if __name__ == "__main__":

    mapFilename='NB-core_v4.map'
    pedFilename='iris.ped'
    pedFilename='IRIS_313-15909.ped'
    genomeFilename='../all.chrs.con'
    gffFilename='../nipponbare-cds-chr.gff3'


    #load genome
    print 'Loading reference genome...'
    genome=[0] ## have first element 0 to shift index
    nchromosomes = 0
    for chrom in SeqIO.parse(genomeFilename, 'fasta', alphabet=generic_dna) :
        genome.append(chrom.seq)
        nchromosomes +=1

    #load annotation
    print 'Loading annotation...'
    # integers, chromosome no, start, end
    annotation=parseGff(gffFilename)
    
    #load mapfile
    print 'Loading map file...'
    mapfile = loadtxt (mapFilename, dtype=int, usecols=(0,3))
    nsnps = mapfile.shape[0]

    #load pedfile
    #pedfile = loadtxt(pedFilename, dtype=str)
    #nsamples = pedfile.shape[0]

    # open pedfile 
    pedfile = open(pedFilename, 'r')
    print 'Parsing ped file...'
    
    for line in pedfile :
        #remove newline and split into columns
        pedline = line[:-1].split(' ')
        
        cultivar, supressed, amplified, stats = pedlineProcess(genome, mapfile, pedline)
        print cultivar, stats 

        newgenome = mutate(genome, mapfile, supressed)

        genes = splitGenome(newgenome, annotation)

        #write new genome to fasta file
        newfasta = open(cultivar+'.fasta', 'w')
        SeqIO.write(genes, newfasta, 'fasta')
        newfasta.close()

    pedfile.close()



