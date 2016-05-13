#functions to calculate distance btw two CDS

import Bio.SubsMat.MatrixInfo as subs
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna, generic_protein
from itertools import izip
from numpy import *
from os.path import isfile
from collections import Counter


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

    #assert len(newseq) %3 == 0, 'cds not of length 3n'
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



# fast hamming distance calculator for any iterable

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


# fast matrix entropy calculator
# returns entropy in bits
def matrix_entropy(matrix) :
    counter = Counter() 
    for i in matrix.flat :
        counter[i] +=1 

    N = float(sum(counter.values()))
    log2N = log2(N)
    
    return -sum([c * (log2(c)-log2N) for c in  counter.itervalues()  ]) / N


blosum = complete_scoring(subs.blosum60)
# calculates a distance matrix and associated metadata
# saves all of them to a npz file

def calc_distmat(gene) :
    gname = gene['Name']
    print gname,
    
    # check if file already exists
    if isfile(gname+'.npz') :
        print 'already exists, skipping it!'
        return None

    # get the full reference sequence
    if gene['strand'] =='+' :   reference = gene['mRNA'].values()[0]['seq-full'].translate()
    else :                      reference = gene['mRNA'].values()[0]['seq-full'].reverse_complement().translate()
    # get the dynamic component of the reference sequence
    reference_dynamic = get_dynamic_for_cultivar(gene, zeros(snps.shape[1],dtype='|S1'))

    #get the score for the static subset of the sequence
    fullscore = score_by_matrix(reference,reference,blosum)[0]
    dynscore = score_by_matrix(reference_dynamic,reference_dynamic,blosum)[0]
    staticScore = fullscore-dynscore
 #   print ' ', fullscore, dynscore, staticScore
    
    nsnps = len(gene['mRNA'].values()[0]['SNPpos'])

    # get the variant proteins for each cultivar
    proteinlist = [ get_dynamic_for_cultivar(gene,pedline) for pedline in snps ]
    # add the reference too
    proteinlist.append(reference_dynamic)
    # go full numpy
    proteinlist=array(proteinlist)
    ncultivars, protsize = proteinlist.shape


    # track if there is a premature termination codon anywhere
    hasPTC = zeros(ncultivars,bool)

    # calculate distances (scores) btw each cultivar
    distmat = zeros ([ncultivars, ncultivars], dtype=int16)

    for i,p1 in enumerate(proteinlist) :
        for j,p2 in enumerate(proteinlist[i:]) :
            dist, hasPTC_this = score_by_matrix(p1,p2,blosum)
            distmat[i][i+j]=dist
            distmat[i+j][i]=dist
        
        hasPTC [i] = hasPTC_this

    # calculate entropy of distance matrix
    s = matrix_entropy(distmat)
    
    # save the distance matrix if it is not trivial
    if s >1e-10 :
        savez_compressed(gname, seq=reference, static=staticScore ,dist=distmat, name=gname, note=gene['Note'], size=protsize, s=s, hasPTC=hasPTC, nsnps=nsnps  )
        print 'has {:.2f} bits of entropy'.format(s)
    else :
        print 'has no entropy, not saving it.'

    # gene name, protein length, entropy, normalized entropy
    return gname, s


if __name__ == '__main__' : 
    from sys import argv
    from glob import iglob

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


    # do not re-do if files exist
    done = [fname.split('.')[0] for fname in iglob('*.npz')]

    todo = []
    for gene in rice.genelist :
        if gene not in done :
            todo.append(gene) 

    print '# going to work on {} remaining genes now.'.format(len(todo))

    from joblib import Parallel, delayed
    log=Parallel(n_jobs=-1,batch_size=100,verbose=10) (delayed(calc_distmat)(rice[gene]) for gene in todo )
#    log=[calc_distmat(gene) for gene in rice ]
    savetxt('complete.log',log,fmt='%s')

        

