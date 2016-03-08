#! /usr/bin/env python2

from sys import argv
from numpy import *
import scipy.cluster.hierarchy as hier
from glob import iglob
from sklearn.metrics import normalized_mutual_info_score
from scipy.spatial.distance import pdist


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

    clusterlist = list()
    clusterlist_noptc = list()

    for fname in iglob('*.npz'):
        data = load(fname)
        hasPTC = data['hasPTC'].any() 
        if hasPTC :
           fullmat = integrate_ptc(data)
        else :
           fullmat = data['dist']

        Z = hier.linkage(fullmat, method='average')

        leaves = hier.leaves_list(Z)

#        newmat=matrix[leaves,:]
#        newmat=newmat[:,leaves]
        
        label = hier.cut_tree(Z, n_clusters=nclusters).flatten()
        
        clusterlist.append(label)
        if not hasPTC : 
            clusterlist_noptc.append(label)

    distmat = pdist(clusterlist, normalized_mutual_info_score)
    savetxt('distmat.dat',distmat)
    distmat_noptc = pdist(clusterlist_noptc, normalized_mutual_info_score)
    savetxt('distmat_noptc.dat',distmat)



