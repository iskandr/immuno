import numpy as np 
import data 
import amino_acid

X,Y = data.load_dataset()

import sklearn
import sklearn.cross_validation 
import sklearn.ensemble

n_classifiers = 5000

clf = sklearn.ensemble.RandomForestClassifier(n_classifiers)

#print "Amino acids"
#print np.mean(sklearn.cross_validation.cross_val_score(clf, X, Y, cv = 10))

fns = [amino_acid.hydropathy, 
       amino_acid.volume, 
       amino_acid.pK_side_chain,
       amino_acid.polarity, 
       amino_acid.prct_exposed_residues,
       amino_acid.hydrophilicity, 
       amino_acid.accessible_surface_area,
       amino_acid.local_flexibility,
       amino_acid.accessible_surface_area_folded,
       amino_acid.refractivity
       ]

#print "All features, all positions"
#X2 = data.transform(X, fns)
#print np.mean(sklearn.cross_validation.cross_val_score(clf, X2, Y, cv = 10))

print "Pairwise ratios"
X2 = data.transform(X, fns, pairwise_ratios = True)
print X2.shape
print np.mean(sklearn.cross_validation.cross_val_score(clf, X2, Y, cv = 10))

#print "Mean per feature"
#X3 = data.transform(X, fns, mean = True)
#print np.mean(sklearn.cross_validation.cross_val_score(clf, X3, Y, cv = 10))

#print "Positions 4,6,8,9"
#X4 = data.transform(X, fns, positions = (4,6,8,9))
#print np.mean(sklearn.cross_validation.cross_val_score(clf, X4, Y, cv = 10))


