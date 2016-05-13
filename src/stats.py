#!/usr/bin/env python2
from parsers import fasta_reference,plink_open,ped_iterator
import gzip



def ped_stats(pedfile, reference=None, verbose=True) :
    
    stats = {}
    firstrow=True

    print '# Parsing PED file :', pedfile.name
    
    for name, seq in ped_iterator(pedfile) :

        if firstrow : 
            print '# nsnps : {:d}'.format(len(seq))
            print '# name\tchanging\thomo\thetero\tmissing'
            firstrow = False
        
        homo   = 0
        hetero = 0
        missing = 0 
        changing = 0

        for i in range(len(seq)) :
            c = seq[i]
            if   c == '0' :  missing += 1
            elif c == ' ' :  hetero  += 1
            elif c in ('A','T','C','G') :
                homo += 1 
                if reference is not None :
                    if reference[i] != c :
                        changing +=1

        if verbose :
            print '{}\t{:d}\t{:d}\t{:d}\t{:d}'.format(name,changing, homo, hetero, missing)

        stats[name]={'homo':homo,'hetero':hetero,'missing':missing, 'changing':changing}

    return stats

if __name__ == '__main__' :

    import argparse

    parser = argparse.ArgumentParser(\
            description='Count changing mutations in a snp set with respect to a reference genome.') 

    parser.add_argument('fastafname', help="fasta file to build the reference")
    parser.add_argument('plinkfname', help="plink file basename for the database")

    args=parser.parse_args()

    mapfile, pedfile = plink_open(args.plinkfname)
    reference = fasta_reference(args.fastafname, mapfile)

    stats = ped_stats(pedfile, reference)

