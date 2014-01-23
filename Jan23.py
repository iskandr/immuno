import numpy as np 
import data 
import amino_acid
import iedb

X,Y = iedb.load_dataset()

import sklearn
import sklearn.cross_validation 
import sklearn.ensemble

n_classifiers = 100

clf = sklearn.ensemble.RandomForestClassifier(n_classifiers)

print "Amino acid histogram vectors"
print np.mean(sklearn.cross_validation.cross_val_score(clf, X, Y, cv = 10))

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
