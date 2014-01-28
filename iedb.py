import numpy as np
import pandas as pd

def load_csv(filename = 'tcell_compact.csv', 
             assay_group=None, 
             unique_sequences = True, 
             filter_noisy_labels = True,
             human = True, 
             hla_type1 = True):
  df = pd.read_csv(filename, skipinitialspace=True)
  mhc = df['MHC Allele Name']

  # 
  # Match known alleles such as 'HLA-A*02:01', 
  # broader groupings such as 'HLA-A2'
  # and unknown alleles of the MHC-1 listed either as 
  #  'HLA-Class I,allele undetermined'
  #  or
  #  'Class I,allele undetermined'
  class_1_mhc_mask = mhc.str.contains('Class I,|HLA-[A-C]([0-9]|\*)', na=False).astype('bool')
  
  print "Class I MHC Entries", class_1_mhc_mask.sum()
  
  # just in case any of the results were from mice or other species, 
  # restrict to humans
  human_mask = df['Host Organism Name'].str.startswith('Homo sapiens', na=False).astype('bool')
  
  print "Human entries", human_mask.sum()
  print "Human Class I MHCs", (human_mask & class_1_mhc_mask).sum()
  
  null_epitope_seq = df['Epitope Linear Sequence'].isnull()
  print "Dropping %d null sequences" % null_epitope_seq.sum()
  # if have rare or unknown amino acids, drop the sequence
  bad_epitope_seq = df['Epitope Linear Sequence'].str.contains('u|x|U|X').astype('bool')
  print "Dropping %d bad sequences" % bad_epitope_seq.sum()
  has_epitope_seq = ~(bad_epitope_seq | null_epitope_seq)
  
  mask = has_epitope_seq
  if human:
    mask &= human_mask
  if hla_type1:
    mask &= class_1_mhc_mask 
  if assay_group:
    mask &= df['Assay Group'] == assay_group
  
  print "Filtered sequences epitope sequences", mask.sum()
  
  df = df[mask]
  
  imm_mask = df['Qualitative Measure'].str.startswith('Positive').astype('bool')
  
  imm = df['Epitope Linear Sequence'][imm_mask]
  print "# immunogenic sequences", len(imm)
  print "sequence length"
  #print imm.str.len().describe()
  
  non_mask = df['Qualitative Measure'] == 'Negative'
  non = df['Epitope Linear Sequence'][non_mask]
  print "# non-immunogenic sequences", len(non)
  #print "sequence length", non.str.len().describe()
  
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
import amino_acid
from amino_acid import letter_to_index

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


def peptide_sequences_to_histogram_vectors(peptides, extra_features = fns):
  n = len(peptides)
  X = np.zeros((n, 20 + len(fns))).astype('float')
  for i, peptide in enumerate(peptides):
    for letter in peptide:
      X[i, letter_to_index(letter)] += 1
    X[i, :] /= len(peptide)
  for fn_idx, fn in enumerate(fns):
    X[i, fn_idx] = np.mean([fn(aa) for aa in peptide])
  print X
  return X

"""
def load_dataset_old(filename = 'tcell_compact.csv', 
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
"""

from sklearn.feature_extraction.text import CountVectorizer
from sklearn.preprocessing import normalize
def load_dataset(filename = 'tcell_compact.csv',
                 assay_group=None,
                 unique_sequences = True,
                 filter_noisy_labels = True,
                 human = True, 
                 hla_type1 = True,
                 max_ngram = 1, 
		         normalize_row = True, 
		         reduced_alphabet = None):
  imm, non = load_csv(filename, 
     assay_group, 
     unique_sequences, 
     filter_noisy_labels, 
     human, 
     hla_type1)
  total = list(imm) + list(non)
  
  if reduced_alphabet is None:
    preprocessor = None
  else:
    def preprocessor(s):
      return ''.join([chr(48 + reduced_alphabet[char]) for char in s])
  
  c = CountVectorizer(analyzer='char', 
                      ngram_range=(1,max_ngram),
                      dtype=np.float, 
                      preprocessor = preprocessor)
  
  
  # returns a sparse matrix 
  X = c.fit_transform(total).todense()
  print "Alphabet", c.get_feature_names()
  if normalize_row:
    X = normalize(X, norm='l1')
  Y = np.ones(len(total), dtype='bool')
  Y[len(imm):] = 0
  return X, Y
  
