from numpy import *

# this file has a bunch of sanity checks. 


# this function checks whether all snps mapped to all CDS intervals are actually in the intervals. 

def snp_map_check(mappedGenome, mapfname) :
    print 'loading file : ', mapfname
    mapraw = loadtxt (mapfname, dtype=int, usecols=(0,3))

    print 'checking snp mapping...'
    count = 0
    for gene in mappedGenome : 
        for mrna in gene['mRNA'].itervalues() :
            for i in range(len(mrna['CDS'])) :
                snplist=mrna['SNP'][i]
                interval = mrna['CDS'][i]
                for snp in snplist :
                    count +=1
                    if mapraw[snp][1] < interval[0] or mapraw[snp][1]>interval[1] :
                        print 'shits all fucked up yo!' 
                        print gene
                        return False

    print '{} mapped SNPs checked against {}'.format(count, mapfname)
    print 'SNP mapping looks sane' 
    return True

# TODO
# implement start-stop codon and 3n length check for sequences

