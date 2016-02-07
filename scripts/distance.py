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




def protein(gene, pedline=None) :
       #take 1st product for simplicity
        product = gene['mRNA'].values()[0]
        
        #return translated product id there is no pedline
        if pedline == None  :
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



from itertools import izip
def hamdist(str1, str2):       
    assert len(str1) == len(str2) 
    diffs = 0
    for ch1, ch2 in izip(str1, str2):
        if ch1 != ch2:
            diffs += 1
    return diffs

#fast entropy calculator
# assuming positive integers and symmetrical matrix
from numpy import *
def matentropy(mat) :

    assert mat.dtype == int
    
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
    from genome import genome
    rice = genome('../data/MSU7/all.chrs.gff3', '../data/MSU7/all.chrs.con', '../data/filtered.map.gz')

    from parsers import PED_parse_homo
    cname, snps = PED_parse_homo('../data/filtered.ped.gz', 100 )
    from pylab import zeros,sum,matshow, show

    entlist=[]
    for gene in rice :

        proteins = [ protein(gene,pedline) for pedline in snps ]

        distmat = zeros ([len(proteins), len(proteins)], dtype=int)

        for i,p1 in enumerate(proteins) :
            for j,p2 in enumerate(proteins[i+1:]) :
                #print i,i+j
                d=hamdist(p1,p2)
                distmat[i][i+j+1]=d

        s = matentropy(distmat)
        print s
        entlist.append[s]
       # if sum(distmat) >len(proteins[0]) : 
         #   matshow(distmat)
          #  show()
