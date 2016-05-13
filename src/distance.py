#functions to calculate distance btw two CDS

import Bio.SubsMat.MatrixInfo as subs
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna, generic_protein

# assembles a product from a given list of cds
def cds_assemble(cdslist) :
    cds = Seq('', generic_dna)
    for s in cdslist :
        cds += s
    return cds


#
def protein(gene, pedline=None) :
       #take 1st product for simplicity
        product = gene['mRNA'].values()[0]
        
        #return translated product id there is no pedline
        if not isinstance(pedline, ndarray)     :
            tmp = cds_assemble(product['seq'])

        else : 
            newseq = [cds.tomutable() for cds in product['seq'] ]
            for i in range(len(newseq)) :
                seq = newseq[i]
                snpind = product['SNPind'][i]
                snppos = product['SNPpos'][i]
                for j in range(len(snppos)) :
           #         print 'yes snp!'
                    ind=snpind[j]
                    pos=snppos[j]
                    newchar = pedline[ind]
                 #   print ind, pos, newchar, 'jio'
                    if newchar not in  ('', '0') :
                        seq[pos] = pedline[ind]

            tmp = cds_assemble(newseq)


        if gene['strand'] == '+' :
            return tmp.translate()
        else :
            return tmp.reverse_complement().translate()


# fast hamming distance calculator
# this should do until we introduce a proper scoring matrix

from itertools import izip
def hamdist(str1, str2):       
    assert len(str1) == len(str2) 
    return sum( ch1 != ch2 for ch1, ch2 in izip(str1, str2))

# very fast hamming distance using numpy
def numhamdist(str1,str2) :
    assert len(str1) == len(str2) 
    assert str1.dtype == '|S1' 
    return sum(str1 != str2)


#fast entropy calculator
# assuming positive integers and symmetrical matrix
from numpy import *
def matentropy(mat) :

    assert mat.dtype == dtype('uint16')
    
    counts = zeros(mat.max()+1)
    for i,line in enumerate(mat) : 
        for d in line[i+1:] :
            counts[d]+=1

    n = len(mat)
    N = (n*(n-1))/2
    logN = log2(N)
    
    entropy = 0.
    for c in counts :
        if c : entropy += c*(log2(c)-logN)

    return -entropy/N
    


if __name__ == '__main__' : 
    from sys import argv

    gff3fname = argv[1]
    fastafname = argv[2]
    mapfname = argv[3]
    pedfname = argv[4]
    
    from genome import genome
    rice = genome(gff3fname, fastafname, mapfname)

    from parsers import ped_parser_homo
    if len(argv)>5 :
        cname, snps = ped_parser_homo(pedfname, int(argv[5]) )
    else : 
        cname, snps = ped_parser_homo(pedfname)
    ncultivars = len(cname)+1

    from numpy import *
    from subprocess import Popen

    for gene in rice :
        # get the variant proteins for each cultivar
        proteinlist = [ array(protein(gene,pedline)) for pedline in snps ]
        # add the reference too
        proteinlist.append(array(protein(gene)) )

        protsize = len(proteinlist[0])
        gname = gene['Name']

        # calculate hamming distances btw each cultivar
        distmat = zeros ([ncultivars, ncultivars], dtype=uint16)

        for i,p1 in enumerate(proteinlist) :
            for j,p2 in enumerate(proteinlist[i+1:]) :
                #d=hamdist(p1,p2)
                # very fast hamming distance using numpy (inline)
                d = sum (p1 != p2) 
                distmat[i][i+j+1]=d
                distmat[i+j+1][i]=d

        # calculate entropy of distance matrix
        s = matentropy(distmat)
        # gene name, protein length, entropy, normalized entropy
        print gname, protsize, s, s/log2(protsize)
        
        # save the distance matrix it is not trivial
        if s >1e-10 :
            save(gname,distmat)
            Popen(['gzip', gname+'.npy'])



