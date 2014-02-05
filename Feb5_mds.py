import scipy 
import scipy.sparse
import numpy as np 

import iedb 
import amino_acid 
from amino_acid import peptide_to_indices

ASSAY = 'cytotoxicity'
LENGTH = 9

imm, non = iedb.load_csv(peptide_length = LENGTH, only_hla_a2 = False, assay_group = ASSAY)
imm = list(imm)
non = list(non)

peptides = imm + non
labels = [True] * len(imm) + [False] * len(non)

X = np.array([peptide_to_indices(p) for p in peptides])
Y = np.array(labels)

n = len(labels)
print "Computing distances"

import sklearn.metrics
import sklearn.metrics.pairwise
D = (sklearn.metrics.pairwise.pairwise_distances(X, metric='hamming') * LENGTH).astype(int)


D = 2 ** (D-1)

print "MDS"
import sklearn.manifold
mds = sklearn.manifold.MDS(dissimilarity="precomputed")

X_2d = mds.fit_transform(D)
print "New dims", X_2d.shape

X_imm_2d = X_2d[:len(imm), :]
X_non_2d = X_2d[len(imm):, :]

import pylab

pylab.scatter(X_imm_2d[:, 0], X_imm_2d[:, 1], c='r')
pylab.scatter(X_non_2d[:, 0], X_non_2d[:, 1], c='g')
pylab.legend(('IMM', 'NON'), loc='best')
pylab.show()