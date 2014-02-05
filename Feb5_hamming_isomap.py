import scipy 
import scipy.sparse
import numpy as np 

import iedb 
import amino_acid 
from amino_acid import peptide_to_indices

CUTOFF = 3
SPARSE = False
ASSAY = None #'cytotoxicity'
LENGTH = 9

imm, non = iedb.load_csv(peptide_length = LENGTH, only_hla_a2 = False, assay_group = ASSAY)
imm = list(imm)
non = list(non)

peptides = imm + non
labels = [True] * len(imm) + [False] * len(non)

X = np.array([peptide_to_indices(p) for p in peptides])
Y = np.array(labels)

n = len(labels)

import sklearn.metrics
import sklearn.metrics.pairwise
D = sklearn.metrics.pairwise.pairwise_distances(X, metric='hamming')
D = np.round(D*LENGTH).astype('int')

print "Distances"  
print D
print (D > 0).sum(axis=1)
print (D > 0).sum()

D_neighbors = D.copy()
D_neighbors[D > CUTOFF] = 0

print "Distances (neighbors)"  
print D_neighbors
D_neighors = scipy.sparse.csr_matrix(D_neighbors)

print (D_neighbors > 0).sum(axis=1)
print (D_neighbors > 0).sum()
print ((D_neighbors > 0).sum(axis=1) == 0).sum()

import sklearn.utils
import sklearn.utils.graph_shortest_path
print "Manifold distances"
D_fw = sklearn.utils.graph_shortest_path.graph_shortest_path(D_neighbors, directed=False, method='FW')
# print D_fw 
print np.sum(D_fw > 0, axis=1)
print np.sum(D_fw > 0)
print (np.sum(D_fw > 0, axis=1) == 0).sum()

import sklearn.manifold
mds = sklearn.manifold.MDS(dissimilarity="precomputed")

X_2d = mds.fit_transform(D_fw)
print "New dims", X_2d.shape

X_imm_2d = X_2d[:len(imm), :]
X_non_2d = X_2d[len(imm):, :]

import pylab

pylab.scatter(X_imm_2d[:, 0], X_imm_2d[:, 1], c='r')
pylab.scatter(X_non_2d[:, 0], X_non_2d[:, 1], c='g')
pylab.legend(('IMM', 'NON'), loc='best')
pylab.show()