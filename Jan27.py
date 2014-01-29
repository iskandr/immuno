import numpy as np 
import data 
import amino_acid
import iedb

X,Y = iedb.load_dataset(assay_group = 'cytotoxicity')

import sklearn
import sklearn.cross_validation 
import sklearn.ensemble
import sklearn.linear_model


print "Amino acid histogram vectors w/ Logistic Regression"
lr = sklearn.linear_model.LogisticRegression()
print "LR Accuracy", np.mean(sklearn.cross_validation.cross_val_score(lr, X, Y, cv = 10))

lr.fit(X,Y)
print "LR coefs", lr.coef_


n_classifiers = 100

rf = sklearn.ensemble.RandomForestClassifier(n_classifiers)

print "Amino acid histogram vectors w/ RF"
print "RF Accuracy", np.mean(sklearn.cross_validation.cross_val_score(rf, X, Y, cv = 10))

rf.fit(X,Y)
print "Features", rf.feature_importances_

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
