_short_names = ["Ala", "Arg", "Asn", "Asp", "Cys", "Glu", "Gln", "Gly", "His", 
	       "Ile", "Leu", "Lys", "Met", "Phe", "Pro", "Ser", "Thr", "Trp", 
	       "Tyr", "Val"]

_letters = ["A", "R", "N", "D", "C", "E", "Q", "G", "H", "I", "L", "K", 
		"M", "F", "P", "S", "T", "W", "Y", "V"]

def letter_to_index(x):
  """
  Convert from an amino acid's letter code to its position index
  """
  assert len(x) == 1
  x = x.upper()
  assert x in _letters
  return _letters.index(x)

 
def peptide_to_indices(xs):
  return [letter_to_index(x) for x in xs]

def letter_to_short_name(x):
  return _short_names[letter_to_index(x)]

def peptide_to_short_names(xs):
  return [letter_to_short_name(x) for x in xs]
 


 
