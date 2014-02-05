
import numpy as np
import imma
from sklearn.cross_validation import cross_val_score
from sklearn.base import BaseEstimator
from sklearn.metrics import roc_auc_score
from collections import defaultdict


"""
Properties of MHC Class I Presented Peptides That
Enhance Immunogenicity

Jorg J. A. Calis1
*, Matt Maybeno2
, Jason A. Greenbaum2
, Daniela Weiskopf2
, Aruna D. De Silva2,3,
Alessandro Sette2
"""

# Table 2 
kl_div = [ 0, 0, 0.10, 0.31, 0.30, 0.29, 0.26, 0.18, 0]

class ImmunoEnhanceModel(BaseEstimator):
  
  def __init__(self):
    self._aa_counts = [defaultdict(int), defaultdict(int)]
  
  def fit(self, X, Y):
    for x, y in zip(X,Y):
      for aa in x:
        self._aa_counts[y][aa] += 1.0


  def predict_proba(self, X):
    return np.array( [ [0, self._predict_instance(x)] for x in X])

  def _predict_instance(self, pep):
    return np.array(np.sum ( [kl_div[i] * 
                                  np.log(self._aa_counts[1][aa] / self._aa_counts[0][aa]
                               )
                    for (i, aa) in enumerate(pep)]))

if __name__ == '__main__':
  model = ImmunoEnhanceModel()
  imm, non = imma.load_csv()
  total = list(imm) + list(non)
  Y = np.ones(len(total), dtype='int')
  Y[len(imm):] = 0
  model.fit(total, Y)
  print model.predict_proba( [ "SLFNTVATL"] )
  #print roc_auc_score(Y, model.predict_proba(total)[:, 1])
  #print np.mean(cross_val_score(model, total, Y, cv=3))
