import numpy as np


import sklearn
import sklearn.cross_validation 
import sklearn.ensemble
import sklearn.linear_model

import data 
import amino_acid
import iedb
import reduced_alphabet 
import eval_dataset


# first try all species
print "All species"
X_all, Y_all = iedb.load_dataset(
                 filter_noisy_labels = False,
                 human = False, 
                 hla_type1 = False,
                 exclude_hla_a2 = False, 
                 only_hla_a2 = False)
eval_dataset.eval_cv(X_all, Y_all)


# also try filtering out the contradictory entries
print
print "All species filtered"
X_all_filter, Y_all_filter = iedb.load_dataset(
                 filter_noisy_labels = True,
                 human = False, 
                 hla_type1 = False,
                 exclude_hla_a2 = False, 
                 only_hla_a2 = False)
eval_dataset.eval_cv(X_all_filter, Y_all_filter)

print "Human"
X_human, Y_human = iedb.load_dataset(
                 filter_noisy_labels = False,
                 human = True, 
                 hla_type1 = False,
                 exclude_hla_a2 = False, 
                 only_hla_a2 = False)
eval_dataset.eval_cv(X_human, Y_human)

print "Human filtered"
X_human_filter, Y_human_filter = iedb.load_dataset(
                 filter_noisy_labels = True,
                 human = True, 
                 hla_type1 = False,
                 exclude_hla_a2 = False, 
                 only_hla_a2 = False)
eval_dataset.eval_cv(X_human_filter, Y_human_filter)

print "Human MHC1"
X_human_mhc1, Y_human_mhc1 = iedb.load_dataset(
                 filter_noisy_labels = False,
                 human = True, 
                 hla_type1 = True,
                 exclude_hla_a2 = False, 
                 only_hla_a2 = False)
eval_dataset.eval_cv(X_human_mhc1, Y_human_mhc1)

print "Human MHC1 filtered"
X_human_mhc1_filter, Y_human_mhc1_filter = iedb.load_dataset(
                 filter_noisy_labels = True,
                 human = True, 
                 hla_type1 = True,
                 exclude_hla_a2 = False, 
                 only_hla_a2 = False)
eval_dataset.eval_cv(X_human_mhc1_filter, Y_human_mhc1_filter)


print "No HLA-A2"
X_no_hla_a2, Y_no_hla_a2 = iedb.load_dataset(
                 filter_noisy_labels = False,
                 human = True, 
                 hla_type1 = True,
                 exclude_hla_a2 = True, 
                 only_hla_a2 = False)
eval_dataset.eval_cv(X_no_hla_a2, Y_no_hla_a2)


print "No HLA-A2 filtered"
X_no_hla_a2_filter, Y_no_hla_a2_filter = iedb.load_dataset(
                 filter_noisy_labels = True,
                 human = True, 
                 hla_type1 = True,
                 exclude_hla_a2 = True, 
                 only_hla_a2 = False)
eval_dataset.eval_cv(X_no_hla_a2_filter, Y_no_hla_a2_filter)


print "Cross-accuracy for HLA-A2 data"
X_hla_a2, Y_hla_a2 = iedb.load_dataset(
                 filter_noisy_labels = True,
                 human = True, 
                 hla_type1 = True,
                 exclude_hla_a2 = False, 
                 only_hla_a2 = True)
eval_dataset.eval_split(X_no_hla_a2, Y_no_hla_a2, X_hla_a2, Y_hla_a2)