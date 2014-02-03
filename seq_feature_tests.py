import numpy as np 
import data 
import amino_acid
import imma, iedb, prop_of_mhc

X,Y = prop_of_mhc.load_dataset()
import sklearn
import sklearn.cross_validation 
import sklearn.ensemble
import sklearn.linear_model

n_classifiers = 100
#clf = sklearn.ensemble.GradientBoostingClassifier(n_estimators = n_classifiers)
clf = sklearn.ensemble.RandomForestClassifier(n_estimators = n_classifiers)
clf.fit(X.todense(), Y)
print "Amino acid histogram vectors"
print np.mean(sklearn.cross_validation.cross_val_score(clf, X.todense(), Y, cv = 2, scoring='roc_auc'))

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
