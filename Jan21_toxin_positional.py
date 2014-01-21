import numpy as np 
import data 
import amino_acid

import sklearn
import sklearn.cross_validation 
import sklearn.ensemble
import sklearn.svm


print "Loading data and transforming to toxin features"
X,Y = data.load_toxin_features(substring_length=3, positional=True)
print X[0]
def run_classifiers(X,Y):
  print "Data shape", X.shape
  for c in [0.0001, 0.001, 0.01, 0.1, 1]:
    svm = sklearn.svm.LinearSVC(C=c)
    print "SVM C =", c
    print np.mean(sklearn.cross_validation.cross_val_score(svm, X, Y, cv = 10))
 
  n_classifiers = 1000
  rf = sklearn.ensemble.RandomForestClassifier(n_classifiers)
  print "Random Forest"
  print np.mean(sklearn.cross_validation.cross_val_score(rf, X, Y, cv = 10))

"""
Toxin features alone seem to do terribly
"""
run_classifiers(X,Y)

