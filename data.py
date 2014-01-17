import numpy as np 
from amino_acid import peptide_to_indices

def load_lines(filename):
  f = open(filename)
  result = f.read().splitlines()
  f.close()
  return result

def load_9mers(filename):
  lines = load_lines(filename)
  result = [peptide_to_indices(line) for line in lines]
  return np.array(result)

def load_dataset(imm_file = 'IMMA2_imm.txt', non_file='IMMA2_non.txt'):
  imm = load_9mers(imm_file) 
  non = load_9mers(non_file)
  X = np.vstack([imm, non])
  Y = np.ones(len(X), dtype='bool')
  Y[len(imm):] = 0
  return X, Y

def transform(X, fns, positions = None, mean = False):
  X2 = []
  for x in X:
    row = []
    for fn in fns:
      row_entries = [fn(xi) for xi in x]
      if positions:
        row_entries = [xi for i, xi in enumerate(row_entries) if i in positions]
      if mean:
        row.append(np.mean(row_entries))
      else:
        row.extend(row_entries)
    X2.append(row)
  X2 = np.array(X2)
  return X2
  