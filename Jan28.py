import numpy as np


import sklearn
import sklearn.cross_validation 
import sklearn.ensemble
import sklearn.linear_model

import data 
import amino_acid
import iedb
import reduced_alphabet 

def run_classifiers(X,Y):
  lr = sklearn.linear_model.LogisticRegression()
  print "LR Accuracy", np.mean(sklearn.cross_validation.cross_val_score(lr, X, Y, cv = 10))

  lr.fit(X,Y)
  print "LR coefs", lr.coef_


  n_classifiers = 200

  rf = sklearn.ensemble.RandomForestClassifier(n_classifiers)
  print "RF Accuracy", np.mean(sklearn.cross_validation.cross_val_score(rf, X, Y, cv = 10))

  rf.fit(X,Y)
  print "RF Features", rf.feature_importances_


print "4 letter alphabet:"
X4,Y4 = iedb.load_dataset(
  assay_group = 'cytotoxicity', 
  reduced_alphabet = reduced_alphabet.gbmr4,
)
run_classifiers(X4, Y4)

print "---"
print 
print "12 letter alphabet:"
X12,Y12 = iedb.load_dataset(
  assay_group = 'cytotoxicity', 
  reduced_alphabet = reduced_alphabet.sdm12,
)

run_classifiers(X12, Y12)

print "---"
print 
print "17 letter alphabet:"
X17,Y17 = iedb.load_dataset(
  assay_group = 'cytotoxicity', 
  reduced_alphabet = reduced_alphabet.hsdm17,
)

run_classifiers(X17, Y17)


print "---"
print 
print "full alphabet:"
X20,Y20 = iedb.load_dataset(
  assay_group = 'cytotoxicity', 
  reduced_alphabet = None,
)
run_classifiers(X20, Y20)



print "4 letter alphabet pairs:"
X4,Y4 = iedb.load_dataset(
  assay_group = 'cytotoxicity', 
  reduced_alphabet = reduced_alphabet.gbmr4,
  max_ngram = 2,
)
run_classifiers(X4, Y4)

print "---"
print 
print "12 letter alphabet pairs:"
X12,Y12 = iedb.load_dataset(
  assay_group = 'cytotoxicity', 
  reduced_alphabet = reduced_alphabet.sdm12, 
  max_ngram = 2,
)

run_classifiers(X12, Y12)

print "---"
print 
print "17 letter alphabet:"
X17,Y17 = iedb.load_dataset(
  assay_group = 'cytotoxicity', 
  reduced_alphabet = reduced_alphabet.hsdm17,
  max_ngram = 2,
)

run_classifiers(X17, Y17)


print "---"
print 
print "full alphabet:"
X20,Y20 = iedb.load_dataset(
  assay_group = 'cytotoxicity', 
  reduced_alphabet = None,
  max_ngram = 2,
)
run_classifiers(X20, Y20)


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
