import gzip
from numpy import *


# get the line for a given cultivar

def ped_find_cultivar( pedfname , cultivar ) :
    f = gzip.open(pedfname,'r')
    
    for line in f :
        if line[:len(cultivar)] == cultivar :
            return line.strip()


def map_find_loci(mapfname, sid, begin, end) :

        mapfile = loadtxt (mapfname, dtype=int, usecols=(0,3))
