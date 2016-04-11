#!/usr/bin/env python2
# coding: utf-8

from numpy import *
import gzip 
import re
from parsers import ped_iterator,strip_vcf,map_parser
from distance import hamdist
from itertools import izip

#find intersection of multidimensional arrays

def multidim_intersect(arr1, arr2):
    arr1_view = arr1.view([('',arr1.dtype)]*arr1.shape[1])
    arr2_view = arr2.view([('',arr2.dtype)]*arr2.shape[1])
    intersected = intersect1d(arr1_view, arr2_view)
    return intersected.view(arr1.dtype).reshape(-1, arr1.shape[1])

def loci_intersect(poslist1, poslist2) :
    return sorted(set(poslist1).intersection(set(poslist2)))

# finds the indices of small array in big 

def find_indices(small, big) :
    index=[]
    j=0
    i=0

    while i <len(small)  :
        while (big[j] < small[i]).any() and j< len(big)-1 :
            j+=1
    
        if (small[i] == big[j]).all():
            index.append(j)
        i += 1
    return index




if __name__ == '__main__' : 

    import argparse,os

    parser = argparse.ArgumentParser() 
    parser.add_argument('plinkfname', help="plink file basename for the database")
    parser.add_argument('vcffname', help="vcf file for the unkown cultivar")

    args=parser.parse_args()
            
    # strip vcf file to get homozygous snps only
    # positions are (chromosome, index) 
    unknownPos, unknownBase = strip_vcf(args.vcffname) 
    unknownName = os.path.splitext(args.vcffname)[0]

    # load plink MAP file 
    # this is the index for SNPs in the PED file
#    mapfile = loadtxt (args.plinkfname+'.map.gz', dtype=int, usecols=(0,3))
    plinkPos = map_parser(args.plinkfname+'.map.gz')

    # get intersection of SNPs in MAP and VCF files
    print 'finding intersection...'
    intersect = multidim_intersect(array(unknownPos),array(plinkPos))
#    intersect = loci_intersect(unknownPos, plinkPos)
    mapindex = find_indices(intersect, plinkPos)
    unknownindex = find_indices(intersect, unknownPos)

    #create reference sequence for the unknown cultivar
    unknownSeq=''
    for i in unknownindex :
        unknownSeq+=unknownBase[i]

    #unknownSeq=unknownBase[unknownindex]

    
    # compare reference sequence with each cultivar in the PED file
    print "calculating differences..."
    namelist=[]
    distlist=[]
    for cultName, cultSeq in ped_iterator(args.plinkfname+'.ped.gz', mapindex) :
        namelist.append(cultName)
        distlist.append(hamdist(unknownSeq, cultSeq)/float(len(intersect)))

    rank = argsort(distlist)

    # write distance list
    with open(unknownName+'.dist', 'w') as f :
        f.write( '# unknown cultivar provided in {}\n'.format(args.vcffname))
        f.write( '# compared against database in {}\n'.format(args.plinkfname))
        f.write( '# intersect\tunknown\tdatabase\n{}\t{}\t{}\n'.format(len(intersect),len(unknownPos),len(plinkPos)))
        for i in range(3) :
            f.write( '# {}\t{:.3f}\t{} \n'.format(i+1, distlist[rank[i]], namelist[rank[i]]))
        for name, dist in izip(namelist,distlist) :
            f.write('{}\t{:.3f}\n'.format(name, dist) )

    # print some stats
    
