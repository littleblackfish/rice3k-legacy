#!/usr/bin/env python2
# coding: utf-8

from numpy import *
import gzip 
import re

#find intersection of multidimensional arrays

def multidim_intersect(arr1, arr2):
    arr1_view = arr1.view([('',arr1.dtype)]*arr1.shape[1])
    arr2_view = arr2.view([('',arr2.dtype)]*arr2.shape[1])
    intersected = intersect1d(arr1_view, arr2_view)
    return intersected.view(arr1.dtype).reshape(-1, arr1.shape[1])


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

# reads ped file line by line and returns a string of SNPs 
# returns only SNPs from the given index
# returns only homozygous SNPs
# returns ' ' for heretozygous SNP at the given indice

def ped_iterator(pedfname, index) :
    f = gzip.open(pedfname, 'r')
    for line in f :
        snpseq=''
        tmp = line.strip().split(' ')
        for i in index :
            alelle1 = tmp[6+2*i]
            alelle2 = tmp[7+2*i]
            if alelle1 == alelle2 :
                snpseq += alelle1 
            else :
                snpseq += ' '
        yield tmp[0], snpseq


# hamming distance

def hamdist(str1, str2):       
    diffs = 0
    for ch1, ch2 in zip(str1, str2):
        if ch1 != ch2:
            diffs += 1
    return diffs

# much faster hamming distance
#
#from itertools import imap
#def hamdist(str1, str2):
#    assert len(str1) == len(str2)
#    return sum(imap(lambda x,y: x!=y, str1, str2))

# reads vsf file line by line, returns stripped version
# only reads homozygous SNPs that are different from the reference
# replaces
# awk 'length($4) == 1 && length($5) == 1 && $5 != "\." { print $1,$2,$4,$5,$10 }' ${input} |\
#      	sed 's/\:.*//' |\
#	sed 's/\(Osj\)\?[Cc]hr0\?//' |\
#	awk '{if ($5 == "1/1") print $1,$2,$3,$4}' 


def strip_vcf (vcffile) :
    f=open(vcffile, 'r')

    
    pos = [] 
    seq = '' 

    for line in f :
        if line[0] != '#' :
            line=line.strip().split()
            if len(line[3]) == 1 and len(line[4]) == 1 and line[4] != '.' :
                genotype = line[9].split(':')[0]
                if genotype == '1/1' :
                    pos .append ( [int(re.search("\d+",line[0]).group()), int(line[1])] )
                    seq +=  ( line[4] )

    return array(pos), seq 
            



from sys import argv

#fname = 'CX133.lean.homo.dat'
#fname = argv[1]
#unknown = loadtxt(fname, usecols=[0,1], dtype = int)
#unknownBase = loadtxt(fname, usecols=[3], dtype = str)
#reference = loadtxt(fname, usecols=[2], dtype = str)


unknownPos, unknownBase = strip_vcf(argv[1]) 

mapfile = loadtxt ('../../data/NB-core_v4/NB-core_v4.map', dtype=int, usecols=(0,3))


# In[22]:

intersect = multidim_intersect(unknownPos,mapfile)
mapindex = find_indices(intersect, mapfile)
unknownindex = find_indices(intersect, unknownPos)


# In[23]:

unknownSeq=''
for i in unknownindex :
    unknownSeq+=unknownBase[i]


# In[31]:

name=[]
dist=[]
for cultivar in ped_iterator('../../data/NB-core_v4/NB-core_v4.ped.gz', mapindex) :
    name.append(cultivar[0])
    dist.append(hamdist(unknownSeq, cultivar[1]))

print argv[1], len(unknownPos), len(mapfile), len(intersect), 1.0-float(min(dist))/float(len(intersect)),  name[argmin(dist)]
