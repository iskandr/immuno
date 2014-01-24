import numpy as np
import pandas as pd

def load_csv(filename = 'tcell_compact.csv', 
             assay_group=None, 
             unique_sequences = True, 
             filter_noisy_labels = True):
  df = pd.read_csv(filename)
  mhc = df['MHC Allele Name']

  # 
  # Match known alleles such as 'HLA-A*02:01', 
  # broader groupings such as 'HLA-A2'
  # and unknown alleles of the MHC-1 listed either as 
  #  'HLA-Class I,allele undetermined'
  #  or
  #  'Class I,allele undetermined'
  class_1_mhc = mhc.str.contains('Class I,|HLA-[A-C]([0-9]|\*)', na=False).astype('bool')
  
  print "Class I MHC Entries", class_1_mhc.sum()
  
  # just in case any of the results were from mice or other species, 
  # restrict to humans
  human = df['Host Organism Name'].str.startswith('Homo sapiens', na=False).astype('bool')
  
  print "Human entries", human.sum()
  print "Human Class I MHCs", (human & class_1_mhc).sum()
  
  null_epitope_seq = df['Epitope Linear Sequence'].isnull()
  print "Dropping %d null sequences" % null_epitope_seq.sum()
  # if have rare or unknown amino acids, drop the sequence
  bad_epitope_seq = df['Epitope Linear Sequence'].str.contains('u|x|U|X').astype('bool')
  print "Dropping %d bad sequences" % bad_epitope_seq.sum()
  has_epitope_seq = ~(bad_epitope_seq | null_epitope_seq)
  
  mask = class_1_mhc & human & has_epitope_seq
  print "Class I MHC humans w/ epitope sequences", mask.sum()
  
  if assay_group:
    mask &= df['Assay Group'] == assay_group
  
  df = df[mask]
  
  imm_mask = df['Qualitative Measure'].str.startswith('Positive').astype('bool')
  
  imm = df['Epitope Linear Sequence'][imm_mask]
  print "# immunogenic sequences", len(imm)
  print "sequence length"
  print imm.str.len().describe()
  
  non_mask = df['Qualitative Measure'] == 'Negative'
  non = df['Epitope Linear Sequence'][non_mask]
  print "# non-immunogenic sequences", len(non)
  print "sequence length", non.str.len().describe()
  
  imm_set = set(imm)
  non_set = set(non)
  noisy_set = imm_set.intersection(non_set)
  print "# unique IMM", len(imm_set)
  print "# unique NON", len(non_set)
  print "# overlap %d (%0.4f)" % (len(noisy_set), float(len(noisy_set)) / len(imm_set))
  
  if unique_sequences:
    if filter_noisy_labels:
      return imm_set.difference(noisy_set), non_set.difference(noisy_set)
    else:
      return imm_set, non_set 
  else:
    if filter_noisy_labels:
      imm_pure = [seq for seq in imm if seq not in noisy_set]
      non_pure = [seq for seq in non if seq not in noisy_set]
      return imm_pure, non_pure
    else:
      return imm, non 
  
import numpy as np 
from amino_acid import letter_to_index
def peptide_sequences_to_histogram_vectors(peptides):
  n = len(peptides)
  X = np.zeros((n, 20)).astype('float')
  for i, peptide in enumerate(peptides):
    for letter in peptide:
      X[i, letter_to_index(letter)] += 1
    X[i, :] /= len(peptide)
  return X
  
def load_dataset(filename = 'tcell_compact.csv', 
                 assay_group=None, 
                 unique_sequences = True, 
                 filter_noisy_labels = True, 
                 balance_classes = True):
  imm, non = load_csv(filename, assay_group, unique_sequences, filter_noisy_labels)
  X_imm = peptide_sequences_to_histogram_vectors(imm)
  X_non = peptide_sequences_to_histogram_vectors(non)
  print "IMM shape", X_imm.shape
  print "NON shape", X_non.shape
  if balance_classes:
    n_imm = len(X_imm)
    n_non = len(X_non)
    if n_imm < n_non:
      n_repeat = int(np.floor(n_non / float(n_imm)))
      X_imms = [X_imm]  * n_repeat 
      X_imm = np.vstack(X_imms)
    else:
      n_repeat = int(np.floor(n_imm / float(n_non)))
      
  X = np.vstack([X_imm, X_non])
  Y = np.ones(len(X), dtype='bool')
  Y[len(X_imm):] = 0
  return X, Y
  