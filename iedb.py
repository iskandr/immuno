import numpy as np
import pandas as pd

def load_csv(filename = 'tcell_compact.csv', 
             assay_group=None, 
             unique_sequences = True, 
             drop_positive_noisy_labels = True,
             drop_negative_noisy_labels = True,
             human = True, 
             hla_type1 = True,
             exclude_hla_a2 = False,
             only_hla_a2 = False,
             nrows = 5000):
  df = pd.read_csv(filename, skipinitialspace=True, nrows = nrows)
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
  bad_epitope_seq = df['Epitope Linear Sequence'].str.contains('u|x|j|b|z|U|X|J|B|Z', na=False).astype('bool')
  print "Dropping %d bad sequences" % bad_epitope_seq.sum()
  has_epitope_seq = ~(bad_epitope_seq | null_epitope_seq)
  
  mask = has_epitope_seq
  if human:
    mask &= human_mask
  if hla_type1:
    mask &= class_1_mhc_mask 
  if assay_group:
    mask &= df['Assay Group'] == assay_group
  
  hla_a2_mask = (mhc == 'HLA-A2') | mhc.str.startswith('HLA-A\*02', na=False)
  print "HLA A-2 count:", hla_a2_mask.sum() 
  if exclude_hla_a2:
    mask &= ~hla_a2_mask
    
  if only_hla_a2:
    mask &= hla_a2_mask
    
  print "Filtered sequences epitope sequences", mask.sum()
  
  df = df[mask]
  
  imm_mask = df['Qualitative Measure'].str.startswith('Positive').astype('bool')
  
  imm = df['Epitope Linear Sequence'][imm_mask]
  print "# immunogenic sequences", len(imm)
  #print "sequence length"
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
    if drop_negative_noisy_labels:
      non_set = non_set.difference(noisy_set)
    if drop_positive_noisy_labels:
      imm_set = imm_set.difference(noisy_set)
    return imm_set, non_set 
  else:
    if drop_negative_noisy_labels:
      non = [seq for seq in non if seq not in noisy_set]
    if drop_positive_noisy_labels:
      imm = [seq for seq in imm if seq not in noisy_set]
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



from sklearn.feature_extraction.text import CountVectorizer
from sklearn.preprocessing import normalize
def load_dataset(filename = 'tcell_compact.csv',
                 assay_group=None,
                 unique_sequences = True,
                 # 'drop' | 'keep' | 'positive' | 'negative'
                 noisy_labels = 'drop', 
                 human = True, 
                 hla_type1 = True,
                 exclude_hla_a2 = False, 
                 only_hla_a2 = False, 
                 max_ngram = 1, 
		         normalize_row = True, 
		         reduced_alphabet = None):
  assert noisy_labels in ('drop', 'keep', 'positive', 'negative'), \
    "Invalid option: %s" % noisy_labels
  drop_positive_noisy_labels = (noisy_labels == 'drop') or (noisy_labels == 'negative')
  drop_negative_noisy_labels = (noisy_labels == 'drop') or (noisy_labels == 'positive')
  
  imm, non = load_csv(filename, 
     assay_group, 
     unique_sequences, 
     drop_positive_noisy_labels,
     drop_negative_noisy_labels, 
     human, 
     hla_type1,
     exclude_hla_a2, 
     only_hla_a2)
  
  print "# IMM", len(imm)
  print "# NON", len(non)
  
  if reduced_alphabet is None:
    preprocessor = None
  else:
    def preprocessor(s):
      return ''.join([chr(48 + reduced_alphabet[char]) for char in s])
  
  c = CountVectorizer(analyzer='char', 
                      ngram_range=(1,max_ngram),
                      dtype=np.float, 
                      preprocessor = preprocessor)
  
  
  total = list(imm) + list(non)
  # returns a sparse matrix 
  X = c.fit_transform(total).todense()
  if reduced_alphabet:
    print "Alphabet", c.get_feature_names()
  if normalize_row:
    X = normalize(X, norm='l1')
  Y = np.ones(len(total), dtype='bool')
  Y[len(imm):] = 0
  print "Dataset size", X.shape
  return X, Y
  
