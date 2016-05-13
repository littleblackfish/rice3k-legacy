import gzip
from numpy import *
from Bio import SeqIO
from Bio.Alphabet import generic_dna, generic_protein
from Bio.SeqRecord import SeqRecord 
from Bio.Seq import Seq
import re
import itertools

def gff3_iterator(fname) :
    f = open(fname, 'r')
    print '# Parsing annotation file :', fname

    for line in f :
        if line[0] == '#' :
            continue
        else : 

            tmp = line.strip().split("\t")

            assert len(tmp) == 9, 'there are not 9 columns in this line, invalid gff3\n{}'.format(line)
            
            # unpack attributes to dictionary
            attrlist = tmp[8].split(";")
            attrdict = dict()

            for attr in attrlist :
                key, value = attr.split("=") 
            #    value= value.split(",")
                attrdict[key]=value

            if 'Note' in attrdict :
                attrdict['Note'] = attrdict['Note'].replace('%2C',',')
                attrdict['Note'] = attrdict['Note'].replace('%20',' ')

            # this is what a generic parser would look like
            # feature = {'seqid':tmp[0], 'source':tmp[1], 'type':tmp[2], 'range':(int(tmp[3]),int(tmp[4])), 'score':tmp[5], 'strand':tmp[6], 'phase':tmp[7], 'attr':attrdict  }

            # but we are making some assumptions

            # seqid is chromosome number, so we convert to int
            seqid = int(filter(str.isdigit,tmp[0])) 
            #sid = int(re.findall('\d+',tmp[0])[0])
            
            # we don't extract source, score and phase because we don't need them
            feature = {'seqid':seqid, 'type':tmp[2], 'range':(int(tmp[3]),int(tmp[4])), 'strand':tmp[6] }
            
            # we flatten by adding attributes to the main dictionary
            feature.update(attrdict)

            yield feature

# parses a gff3 to a dictionary
# hierarchy goes like chromosome - gene - mrna - cds,exon,intron,utr
def gff3_parser(fname) :
    genome=dict()
    for f in gff3_iterator(fname) :

        # we don't use source or score for this application
        #f.pop('source')
        #f.pop('score')
        
        typ = f.pop('type')
        sid = f['seqid']
        fid = f.pop('ID', None)
        
        if typ == 'gene' :
            #f.pop('phase') 

            if sid not in genome :
                genome[sid] = dict()
            
            f['mRNA']=dict()

            genome[sid][fid] =  f
        else : 
            
            parent = f.pop('Parent')

            if typ == 'mRNA' :
                #f.pop('phase')
                f['CDS'] = list()
                f['SNP'] = list()
                #f['exon'] = list()
                #f.pop('attr')
                f.pop('strand')
                f.pop('range')
                genome[sid][parent]['mRNA'][int(fid.split('.')[1])] = f

            else :
                if typ == 'CDS' :
                    gene = parent.split('.')[0]
                    no = int(parent.split('.')[1])

                    genome[sid][gene]['mRNA'][no][typ].append(f['range'])

                if typ == 'exon' : 
                    pass
                if typ == 'five_prime_UTR' :
                    pass
                if typ == 'three_prime_UTR' :
                    pass
    return genome
            

# parse a genome fasta file and return a dictionary of chromosomes
# keys (chromosome no) are integers for convenience 

def fasta_parser(fname) : 
    print '# Parsing fasta file :', fname
    genome = []

    for chrom in SeqIO.parse(fname, 'fasta', alphabet=generic_dna) :
        genome.append(chrom.seq)

    genomedict = {i+1:genome[i] for i in range(len(genome))}

    return genomedict


def fasta_reference(fastafname, mapfile) :
    genome  = fasta_parser(fastafname)
    maplist = map_parser(mapfile)
    
    reference = [genome[m[0]][m[1]-1] for m in maplist]

    return ''.join(reference)


# convenience function to open coupled map/ped files

def plink_open(basename) :
    try :       mapfile = open(basename+'.map') 
    except :    mapfile = gzip.open(basename+'.map.gz')
    
    try :       pedfile = open(basename+'.ped')
    except :    pedfile = gzip.open(basename+'.ped.gz')

    return mapfile, pedfile

def map_parser(mapfile) :
    print '# Parsing MAP file :', mapfile.name
    mapraw = loadtxt (mapfile, dtype=int, usecols=(0,3))

    return sorted([(line[0], line[1]) for line in mapraw])

# parses a MAP file into a dictionary 

def map_dict(mapfile) : 
    print '# Parsing MAP file :', mapfile.name
    mapraw = loadtxt (mapfile, dtype=int, usecols=(0,3))

    mapdict = { i:[] for i in range(1,13)}
    for i in range(len(mapraw)) :
        mapdict[mapraw[i][0]].append([mapraw[i][1], i])

    # make sure it is sorted by position
    for i in range(1,13) :
        mapdict[i]=array(sorted(mapdict[i], key=lambda x:x[0]) ).T
        #mapdict[i]=array(mapdict[i] ).T
    
    return mapdict

# returns SNP indices and positions for a given (closed) interval

def map_find_loci(mapdict, sid, interval) :
    pos = mapdict[sid][0]
    indices = where((pos>=interval[0]) & (pos<=interval[1])) [0]
    lociSNPind = [mapdict[sid][1][i] for i in indices ]
    lociSNPpos = [mapdict[sid][0][i] for i in indices ]
    return lociSNPind, lociSNPpos


# reads ped file line by line and returns a string of SNPs 
# returns only homozygous SNPs (including '0')
# returns ' ' for heretozygous 
# optionally returns only SNPs from the given index

from itertools import izip,islice

# extract homozygous snps from a ped line
# this includes missing (0) calls
# replaces hetero snps with 

def pedline_homo(line) :
    # split line
    line = line.strip().split()
    
    # first column is sample name
    name = line[0]

    # split line into two strands
    # (be careful, data may not be phased)
    strand1 = islice(line, 6, None, 2)
    strand2 = islice(line, 7, None, 2)
    
    seq = zeros((len(line)-6) /2, dtype='|S1')
   
    for i, (alelle1, alelle2) in enumerate(izip(strand1,strand2)) :
        # homozgyous
        if alelle1 == alelle2 : seq[i] = alelle1 
        # heterozygous
        else : seq[i] = ' '

    return name, seq

# iterates over a ped file returning homozygous 

def ped_iterator(pedfile, index=None) :

    for line in pedfile :
        # read line
        # first column is sample name

        name, seq = pedline_homo(line)
        
        # optionally filter for index
        if index is not None :
            seq = seq[index]
        
        yield name, seq

# get the line for a given cultivar
# much faster than iterator because does not split every line
# identical output

def ped_find_cultivar( pedfile , cultivar ) :
    print '# Parsing PED file :', pedfile.name
    
    for line in pedfile :
        if line[:len(cultivar)] == cultivar :
            
            name, seq = pedline_homo(line)

            pedfile.seek(0)
            return name, seq

# parses a ped file 
# puts it in a '|S1' numpy array
# returns a name list and the array
# this is as memory efficient as it gets 
# still requires the whole thing fits into memory though

def ped_parser_homo(pedfile, nrows=3023, index=None) :
    print '# Parsing PED file :', pedfile.name
    
    snps=None
    names = []
    row=0

    for name, seq in ped_iterator(pedfile, index=index) :
        
        # initialize matrix after reading first line
        # this way you know nsnps
        if  snps is None:
            nsnps = (len(seq))
            snps = zeros([nrows,nsnps], dtype='S1')
        
        names.append(name)
        
        snps[row] = seq 
        row +=1
        if row==nrows : break
    
    print '# {} SNPs in {} cultivars.'.format(nsnps, nrows )
    return names,snps


# reads vsf file line by line, returns stripped version
# only reads homozygous SNPs that are different from the reference
# grandfather of this was something like
# awk 'length($4) == 1 && length($5) == 1 && $5 != "\." { print $1,$2,$4,$5,$10 }' ${input} |\
#      	sed 's/\:.*//' |\
#	sed 's/\(Osj\)\?[Cc]hr0\?//' |\
#	awk '{if ($5 == "1/1") print $1,$2,$3,$4}' 
# returns (chromosome, index), sequence


def strip_vcf (vcffile) :
    with open(vcffile, 'r') as f:

        pos = [] 
        seq = [] 

        for line in f :
            if line[0] != '#' :
                line=line.strip().split()
                if len(line[3]) == 1 and len(line[4]) == 1 and line[4] != '.' :
                    genotype = line[9].split(':')[0]
                    if genotype == '1/1' :
                        # extract numerical chromosome no
                        chrno= int(re.search("\d+",line[0]).group())
                        index= int(line[1])
                        pos.append ( (chrno, index) )
                        seq.append ( line[4] )

    return pos, seq
            

# better vcf parser using pyvcf

import vcf

def vcf_homo_snp (vcffile, minQual=50) : 
    reader = vcf.Reader(vcffile) 

    nsamples = len(reader.samples)

    assert nsamples == 1, 'function only implemented for single sample vcf files so far\n{}'.format(vcffile.name) 

    sample = reader.samples[0]

    positions = []
    sequence = ''
    for record in reader : 
        if not record.is_monomorphic and record.is_snp and record.QUAL>minQual :
            g1,g2 = record.genotype(sample)['GT'].split('/')

            if g1==g2 :
                g1=int(g1)
                assert g1 == 1, '{} genotype at {},{}'.format(record.genotype(sample)['GT'], record.CHROM, record.POS)

                chrno =int(re.search("\d+",record.CHROM).group())
                positions.append((chrno, record.POS))
                
                # weird struct here but ok.
                # maybe polish later

                assert len(record.ALT[0].sequence) == 1, 'ALT is weird.'
                
                sequence += record.ALT[0].sequence

    return positions, sequence

