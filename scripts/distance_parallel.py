#functions to calculate distance btw two CDS

import Bio.SubsMat.MatrixInfo as subs
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna, generic_protein
from itertools import izip


# returns the 'dynamic' subset of the translated product
# dynamic meaning those aminoacids with snps on them

def get_dynamic_for_cultivar(gene, pedline) :
   #take 1st mrna for simplicity
    mrna = gene['mRNA'].values()[0]
    
    newseq = mrna['seq-full'].tomutable()

    for snpind,snppos in izip(mrna['SNPind'],mrna['SNPpos']) :
        newchar = pedline[snpind]
        if newchar not in  ('', '0') :
            newseq[snppos] = newchar

    assert len(newseq) %3 == 0, 'cds not of length 3n'
    newseq = newseq.toseq()

    if gene['strand'] == '+' :
        product  = newseq.translate()
        # select 'dynamic' AAs and remove duplicates
        dynamics = sorted(set([i/3 for i in mrna['SNPpos']] ))
        
    else :
        product = newseq.reverse_complement().translate()
        # re-map snp positions for reverse complement 
        flipSNPpos = [len(newseq) - i - 1   for i in mrna['SNPpos']]
        dynamics = sorted(set([i/3 for i in flipSNPpos]))

    return [product[i] for i in dynamics ]



# fast hamming distance calculator
# this should do until we introduce a proper scoring matrix

def hamdist_iter(str1, str2):       
    assert len(str1) == len(str2) 
    return sum( ch1 != ch2 for ch1, ch2 in izip(str1, str2))

# very fast hamming distance using numpy
def hamdist_numpy(str1,str2) :
    assert len(str1) == len(str2) 
    assert str1.dtype == '|S1' 
    return sum(str1 != str2)


# completes a scoring matrix by filling in symmetrical fields 
def complete_scoring (smat) :
    fullmat = dict()
    for aa1,aa2 in smat.keys() :
        fullmat[(aa1,aa2)] = smat[(aa1,aa2)]
        fullmat[(aa2,aa1)] = smat[(aa1,aa2)]
    # add stop codon
    fullmat[('*','*')]=0
    return fullmat

# computes distance using a scoring matrix
def score_by_matrix(seq1, seq2, scoring_matrix) :
    assert len(seq1)==len(seq2), 'comparing sequences of different length'
    score = 0 
    hasPTC = False
    for pair in izip(seq1,seq2) :
        if (pair[0] == '*') ^ (pair[1] == '*') :
            hasPTC=True
            break
        else : 
            score += scoring_matrix[pair] 
    return score, hasPTC


#fast entropy calculator
# assuming positive integers and symmetrical matrix
from numpy import *
def matrix_entropy(matrix) :

    assert matrix.dtype == dtype('int16')
    
    counts = zeros(matrix.max()+1)
    for i,line in enumerate(matrix) : 
        for d in line[i+1:] :
            counts[d]+=1

    n = len(matrix)
    N = (n*(n-1))/2
    logN = log2(N)
    
    entropy = 0.
    for c in counts :
        if c : entropy += c*(log2(c)-logN)

    return -entropy/N

from numpy import *
from subprocess import Popen
from os.path import isfile

blosum = complete_scoring(subs.blosum60)
def calc_distmat(gene) :
    gname = gene['Name']
    print gname,
    
    # check if file already exists
    if isfile(gname+'.npz') :
        print 'already exists, skipping it!'
        return None

    # get the variant proteins for each cultivar
    proteinlist = [ get_dynamic_for_cultivar(gene,pedline) for pedline in snps ]
    # add the reference too
#    proteinlist.append(protein(gene) )
    # go full numpy
    proteinlist=array(proteinlist)
    
    # length of this protein
    protsize = proteinlist.shape[1]

    # track if there is a premature termination codon anywhere
    hasPTC = 0

    # calculate distances (scores) btw each cultivar
    distmat = zeros ([ncultivars, ncultivars], dtype=int16)

    for i,p1 in enumerate(proteinlist) :
        for j,p2 in enumerate(proteinlist[i:]) :
            dist, hasPTC_this = score_by_matrix(p1,p2,blosum)
            distmat[i][i+j]=dist
            distmat[i+j][i]=dist
        
        hasPTC += hasPTC_this

    distmat -= distmat.min()
    # calculate entropy of distance matrix
    s = matrix_entropy(distmat)
    
    # save the distance matrix it is not trivial
    if s >1e-10 :
        savez_compressed(gname, dist=distmat, name=gname, note=gene['Note'], size=protsize, s=s, hasPTC=hasPTC  )
        print 'has {:.2f} bits of entropy and {} PTCs.'.format(s,hasPTC)
    else :
        print 'is trivial (s=0),  it.'

    # gene name, protein length, entropy, normalized entropy
    return gname, protsize, s, s/log2(protsize)

    


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


    from joblib import Parallel, delayed
  #  a= Parallel(n_jobs=-1,verbose=5) (delayed(calc_distmat)(gene) for gene in rice )
    [calc_distmat(gene) for gene in rice ]

        

