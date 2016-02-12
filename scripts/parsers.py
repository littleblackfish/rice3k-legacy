import gzip
from numpy import *
from Bio import SeqIO
from Bio.Alphabet import generic_dna, generic_protein
from Bio.SeqRecord import SeqRecord 
from Bio.Seq import Seq
from parsers import *




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
    print '# Parsing sequence file :', fname
    genome = []

    for chrom in SeqIO.parse(fname, 'fasta', alphabet=generic_dna) :
        genome.append(chrom.seq)

    genomedict = {i+1:genome[i] for i in range(len(genome))}

    return genomedict


# get the line for a given cultivar

def ped_find_cultivar( pedfname , cultivar ) :
    f = gzip.open(pedfname,'r')
    
    for line in f :
        if line[:len(cultivar)] == cultivar :
            return line.strip()

# parses a MAP file into a dictionary 

def map_parser(fname) :
    print '# Parsing MAP file :', fname
    mapdict = { i:[] for i in range(1,13)}
    mapraw = loadtxt (fname, dtype=int, usecols=(0,3))

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

# reads a ped file, takes only homozygous snps 
# returns them in a '1S' numpy array
# along with a name list
# as memory efficient as it gets 

def ped_parser_homo(pedfname, nrows=3023) :
    f = gzip.open(pedfname, 'r') 
    print '# Parsing PED file :', pedfname
    
    ncols = len(f.readline().strip().split(' '))
    nsnps = (ncols-6) / 2
    f.rewind()
    
    snps = zeros([nrows,nsnps], dtype='S1')
    names = []

    row=0
    for line in f :
        line = line.strip().split(' ')
        names.append(line[0])
        for i in range(6, len(line), 2) :
            if line[i] == line[i+1] :
                snps[row,(i-6)/2] = line[i] 
        row +=1
        if row==nrows :
            break
    
    print '# {} SNPs in {} cultivars.'.format(nsnps, nrows )
    return names,snps

def ped_stats(pedfname, nrows=3023) :
    f = gzip.open(pedfname, 'r') 
    print '# Parsing PED file :', pedfname
    
    # figure out number of snps
    ncols = len(f.readline().strip().split(' '))
    nsnps = (ncols-6) / 2
    f.rewind()
    
    stats = {}

    for row,line in enumerate(f) :
        line = line.strip().split(' ')
        name = line[0] #cultivar name
        homoCount   = 0
        heteroCount = 0
        missingCount = 0 
        for i in range(6, len(line), 2) :
            assert line[i] in ('A','T','C','G','0'), 'wtf is this {}'.format(line[i])
            if line[i] == '0' :
                missingCount+=1
            elif line[i] == line[i+1] :
                homoCount +=1
            else :
                heteroCount +=1 
        stats[name]=(homoCount,heteroCount,missingCount)

        if row==nrows :
            break

    return stats



