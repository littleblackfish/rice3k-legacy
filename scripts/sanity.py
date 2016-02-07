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
                snplist=mrna['SNPind'][i]
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

def cds_check(ingenome) :
    count = 0 
    for protein in ingenome.protein_iter() :
        if protein[:3] != 'ATG' : 
            print 'No start codon in CDS {}, got {} instead'.format(count,protein[:3])
            return False
        if len(protein)%3 != 0 :
            print 'CDS {} not of length 3N, it has {} extra,'.format(count, len(protein)%3) ,
            shift = len(protein)%3
            if protein[-3-shift:-shift] not in ('TAG', 'TAA', 'TGA') :                
                print 'No stop codon in CDS {}, got {} instead.'.format(count, protein[-3-len(protein)%3:])
            else : 
                print 'but last proper codon is a stop codon' 

           # return False
        elif protein[-3:] not in ('TAG', 'TAA', 'TGA') :
            print 'No stop codon in CDS {}, got {} instead.'.format(count, protein[-3-len(protein)%3:])
           # return False
        
        count += 1
    
    print '{} cds products checked for start, stop and length'.format(count)
    return True


if __name__ == '__main__' :
    
    gff3fname = '../data/MSU7/all.chrs.gff3'
    fastafname = '../data/MSU7/all.chrs.con'
    mapfname = '../data/NB-core_v4/NB-core_v4.map.gz'
    mapfname = '../data/NB-core_v4/filtered.map.gz'

    from genome import genome
    rice = genome(gff3fname, fastafname, mapfname) 

    snp_map_check(rice, mapfname) 

    cds_check(rice)

