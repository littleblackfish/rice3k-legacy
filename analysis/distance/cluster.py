#! /usr/bin/env python2

from sys import argv
from numpy import *
import scipy.cluster.hierarchy as hier
from glob import iglob,glob
from sklearn.metrics import normalized_mutual_info_score
from scipy.spatial.distance import pdist,squareform


def cluster(matrix) : 
    Z = hier.linkage(matrix, method='average')
    
    leaves = hier.leaves_list(Z)

    newmat=matrix[leaves,:]
    newmat=newmat[:,leaves]

    return leaves, newmat

def integrate_ptc(data) :
    smat = data['dist']
    ptc  = data['hasPTC']
    ptc = (1-ptc) * int(smat.mean())
    fullmat = vstack([smat,ptc])
    return fullmat.T
    


if __name__=='__main__' :

    nclusters=10
    nsamples = 100

    clusterlist = list()
    clusterlist_noptc = list()
    
    filelist = glob('blosum60/*.npz') 
    from random import sample
   

    for fname in sample(filelist, nsamples):
        data = load(fname)
	print fname
        hasPTC = data['hasPTC'].any() 
        if hasPTC :
           fullmat = integrate_ptc(data)
        else :
           fullmat = data['dist']

        Z = hier.linkage(fullmat, method='centroid')

        leaves = hier.leaves_list(Z)
#        newmat=fullmat[leaves,:]
#        newmat=newmat[:,leaves]
        
        label = hier.fcluster(Z, nclusters, criterion='maxclust')
        
        clusterlist.append(label)
        if not hasPTC : 
            clusterlist_noptc.append(label)

    distmat = pdist(clusterlist, normalized_mutual_info_score)
    savetxt('cluster_nmi-{}.dat'.format(nsamples), squareform(distmat))
    distmat_noptc = pdist(clusterlist_noptc, normalized_mutual_info_score)
    savetxt('cluster_nmi_noptc-{}.dat'.format(nsamples),squareform(distmat_noptc))




