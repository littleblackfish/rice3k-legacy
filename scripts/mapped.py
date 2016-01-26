import gzip
from numpy import *


# get the line for a given cultivar

def ped_find_cultivar( pedfname , cultivar ) :
    f = gzip.open(pedfname,'r')
    
    for line in f :
        if line[:len(cultivar)] == cultivar :
            return line.strip()

# parses a MAP file into a dictionary 

def parse_map(mapfile) :
    mapdict = { i:[] for i in range(1,13)}
    mapraw = loadtxt (mapfile, dtype=int, usecols=(0,3))

    for i in range(len(mapraw)) :
        mapdict[mapraw[i][0]].append([mapraw[i][1], i])


    # make sure it is sorted by position
    for i in range(1,13) :
        mapdict[i]=array(sorted(mapdict[i], key=lambda x:x[0]) ).T
        #mapdict[i]=array(mapdict[i] ).T
    
    return mapdict



# returns snp indices for a given (closed) interval
# mapfile = loadtxt (mapfname, dtype=int, usecols=(0,3))

#def map_find_loci(mapdict, sid, interval) :
#    locisnps=list()
#    snplist = mapdict[sid]
#    for snp in snplist :
#        if interval[0] <= snp[0] <= interval[1] :
#            locisnps.append(snp[1])
#
#    return locisnps

def map_find_loci(mapdict, sid, interval) :
    pos = mapdict[sid][0]
    indices = where((pos>=interval[0]) & (pos<=interval[1])) [0]
    locisnps = [mapdict[sid][1][i] for i in indices ]
    return locisnps

# reads a ped file, takes only homozygous snps 
# returns them in a '1S' numpy array
# as memory efficient as it gets 

def homo_ped(pedfname, nrows=3023) :
    f = gzip.open(pedfname, 'r') 
    
    ncols = len(f.readline().strip().split(' '))
    nsnps = (ncols-6) / 2
    f.rewind()
    snp= zeros([nrows,nsnps], dtype='S1')

    row=0
    for line in f :
        line = line.strip().split(' ')
        for i in range(6, len(line), 2) :
            if line[i] == line[i+1] :
                snp[row,(i-6)/2] = line[i] 
        row +=1
        if row==nrows :
            break
    
    print nsnps, nrows 
    return snp


