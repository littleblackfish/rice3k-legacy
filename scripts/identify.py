#!/usr/bin/env python2
# coding: utf-8

from numpy import *
import gzip 
import re
from parsers import ped_iterator,map_parser,vcf_homo_snp,plink_open
from distance import hamdist
from itertools import izip


def pos_intersect(poslist1, poslist2) :
    
    intersection = sorted(set(poslist1).intersection(set(poslist2)))
    
    index1=[]
    index2=[]

    for i,pos in enumerate(poslist1) :
        if pos in intersection :
            index1.append(i)

    for i,pos in enumerate(poslist2) :
        if pos in intersection : 
            index2.append(i)

    return index1, index2, intersection

if __name__ == '__main__' : 

    import argparse,os

    parser = argparse.ArgumentParser() 
    parser.add_argument('plinkfname', help="plink file basename for the database")
    parser.add_argument('vcffname', help="vcf file for the unkown cultivar")
    parser.add_argument('-v', help="print as you calculate",
                    action="store_true")

    args=parser.parse_args()
            
    # strip vcf file to get homozygous snps only
    # positions are (chromosome, index) 
    with open(args.vcffname) as vcffile :
        vcfpos, vcfseq = vcf_homo_snp(vcffile)

    unknownName = os.path.splitext(args.vcffname)[0]

    # load plink MAP file 

    mapfile,pedfile = plink_open(args.plinkfname) 

    # this is the index for SNPs in the PED file
    plinkpos = map_parser(mapfile)

    # get intersection of SNPs in MAP and VCF files
    print '# Finding intersection...'
    vcfindex, mapindex, intersection = pos_intersect(vcfpos, plinkpos)

    #create reference sequence for the unknown cultivar
    vcfintseq=''
    for i in vcfindex :
        vcfintseq += vcfseq[i]

    # compare reference sequence with each cultivar in the PED file
    print "# Calculating hamming distances..."
    namelist=[]
    distlist=[]
    for cultName, cultSeq in ped_iterator(pedfile, index=mapindex) :
        namelist.append(cultName)
        dist = hamdist(vcfintseq, cultSeq)/float(len(intersection))
        distlist.append(dist)
        if args.v :
            print cultName, dist

    rank = argsort(distlist)

    # write distance list
    with open(unknownName+'.dist', 'w') as f :
        f.write( '# unknown cultivar provided in {}\n'.format(args.vcffname))
        f.write( '# compared against database in {}\n'.format(args.plinkfname))
        f.write( '# intersect\tunknown\tdatabase\n{}\t{}\t{}\n'.format(len(intersection),len(vcfpos),len(mappos)))
        for i in range(3) :
            f.write( '# {}\t{:.3f}\t{} \n'.format(i+1, distlist[rank[i]], namelist[rank[i]]))
        for name, dist in izip(namelist,distlist) :
            f.write('{}\t{:.3f}\n'.format(name, dist) )

    # print some stats
    
