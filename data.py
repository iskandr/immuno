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
  return load_9mers(imm_file), load_9mers(non_file)
