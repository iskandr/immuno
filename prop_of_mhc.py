from sklearn.feature_extraction.text import CountVectorizer
from sklearn.preprocessing import normalize
import pandas as pd
import numpy as np

def load_csv(imm_file = 's1.csv'):
  data = pd.read_csv(imm_file)
  if imm_file == 's1.csv':
    imm = data[ (data.Immunogenicity == 'immunogenic') & (data.Species == 'Mus')].Peptide
    non = data[ (data.Immunogenicity == 'non-immunogenic') & (data.Species == 'Mus')].Peptide
  elif imm_file == 's2.csv':
    imm = data[ (data['epitope status'] == 'epitope') & (data.host == 'Homo')].peptide
    non = data[ (data['epitope status'] == 'non-epitope') & (data.host == 'Homo')].peptide
  return imm, non
  

def load_dataset(normalizeRow = False):
  c = CountVectorizer(analyzer='char', ngram_range=(1,2), dtype=np.float)
  imm, non = load_csv()
  total = list(imm) + list(non)
  X = c.fit_transform(total)
  if normalizeRow:
  	X = normalize(X, norm='l1')
  Y = np.ones(len(total), dtype='bool')
  Y[len(imm):] = 0
  return X, Y
