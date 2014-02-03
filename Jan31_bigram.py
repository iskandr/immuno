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

max_ngram = 2
"""
# first try all species
print "All species"
X_all, Y_all = iedb.load_dataset(
                 noisy_labels = 'keep',
                 human = False, 
                 hla_type1 = False,
                 exclude_hla_a2 = False, 
                 only_hla_a2 = False,
                 max_ngram = max_ngram)
eval_dataset.eval_cv(X_all, Y_all)


# also try filtering out the contradictory entries
print
print "---"
print "All species filtered"
X_all_filter, Y_all_filter = iedb.load_dataset(
                 noisy_labels = 'drop',
                 human = False, 
                 hla_type1 = False,
                 exclude_hla_a2 = False, 
                 only_hla_a2 = False,
                 max_ngram = max_ngram)
eval_dataset.eval_cv(X_all_filter, Y_all_filter)


# also try filtering out the contradictory entries
print
print "---"
print "All species noisy = positive"
X_all_positive, Y_all_positive = iedb.load_dataset(
                 noisy_labels = 'positive',
                 human = False, 
                 hla_type1 = False,
                 exclude_hla_a2 = False, 
                 only_hla_a2 = False,
                max_ngram = max_ngram)
eval_dataset.eval_cv(X_all_positive, Y_all_positive)


# also try filtering out the contradictory entries
print
print "---"
print "All species noisy = negative"
X_all_negative, Y_all_negative = iedb.load_dataset(
                 noisy_labels = 'negative',
                 human = False, 
                 hla_type1 = False,
                 exclude_hla_a2 = False, 
                 only_hla_a2 = False,
                 max_ngram = max_ngram)
eval_dataset.eval_cv(X_all_negative, Y_all_negative)

print 
print "---"
print "Human"
X_human, Y_human = iedb.load_dataset(
                 noisy_labels = 'keep',
                 human = True, 
                 hla_type1 = False,
                 exclude_hla_a2 = False, 
                 only_hla_a2 = False,
                 max_ngram = max_ngram)
eval_dataset.eval_cv(X_human, Y_human)


print 
print "---"
print "Human filtered"
X_human_filter, Y_human_filter = iedb.load_dataset(
                 noisy_labels = 'drop',
                 human = True, 
                 hla_type1 = False,
                 exclude_hla_a2 = False, 
                 only_hla_a2 = False,
                 max_ngram = max_ngram)
eval_dataset.eval_cv(X_human_filter, Y_human_filter)


print 
print "---"
print "Human noisy = positive"
X_human_positive, Y_human_positive = iedb.load_dataset(
                 noisy_labels = 'positive',
                 human = True, 
                 hla_type1 = False,
                 exclude_hla_a2 = False, 
                 only_hla_a2 = False,
                 max_ngram = max_ngram)
eval_dataset.eval_cv(X_human_positive, Y_human_positive)


print 
print "---"
print "Human noisy = negative"
X_human_negative, Y_human_negative = iedb.load_dataset(
                 noisy_labels = 'negative',
                 human = True, 
                 hla_type1 = False,
                 exclude_hla_a2 = False, 
                 only_hla_a2 = False,
                 max_ngram = max_ngram)
eval_dataset.eval_cv(X_human_negative, Y_human_negative)
"""
print 
print "---"
print "Human MHC1"
X_human_mhc1, Y_human_mhc1 = iedb.load_dataset(
                 noisy_labels = 'keep',
                 human = True, 
                 hla_type1 = True,
                 exclude_hla_a2 = False, 
                 only_hla_a2 = False,
                 max_ngram = max_ngram)
eval_dataset.eval_cv(X_human_mhc1, Y_human_mhc1)


print 
print "---"
print "Human MHC1 filtered"
X_human_mhc1_filter, Y_human_mhc1_filter = iedb.load_dataset(
                 noisy_labels = 'drop',
                 human = True, 
                 hla_type1 = True,
                 exclude_hla_a2 = False, 
                 only_hla_a2 = False,
                 max_ngram = max_ngram)
eval_dataset.eval_cv(X_human_mhc1_filter, Y_human_mhc1_filter)


print 
print "---"
print "Human MHC1 noisy = positive"
X_human_mhc1_positive, Y_human_mhc1_positive = iedb.load_dataset(
                 noisy_labels = 'positive',
                 human = True, 
                 hla_type1 = True,
                 exclude_hla_a2 = False, 
                 only_hla_a2 = False,
                 max_ngram = max_ngram)
eval_dataset.eval_cv(X_human_mhc1_positive, Y_human_mhc1_positive)

print 
print "---"
print "Human MHC1 noisy = negative"
X_human_mhc1_negative, Y_human_mhc1_negative = iedb.load_dataset(
                 noisy_labels = 'negative',
                 human = True, 
                 hla_type1 = True,
                 exclude_hla_a2 = False, 
                 only_hla_a2 = False,
                 max_ngram = max_ngram)
eval_dataset.eval_cv(X_human_mhc1_positive, Y_human_mhc1_positive)

print 
print "---"
print "No HLA-A2"
X_no_hla_a2, Y_no_hla_a2 = iedb.load_dataset(
                 noisy_labels = 'keep',
                 human = True, 
                 hla_type1 = True,
                 exclude_hla_a2 = True, 
                 only_hla_a2 = False,
                 max_ngram = max_ngram)
eval_dataset.eval_cv(X_no_hla_a2, Y_no_hla_a2)


print 
print "---"
print "No HLA-A2 filtered"
X_no_hla_a2_filter, Y_no_hla_a2_filter = iedb.load_dataset(
                 noisy_labels = 'drop',
                 human = True, 
                 hla_type1 = True,
                 exclude_hla_a2 = True, 
                 only_hla_a2 = False,
                 max_ngram = max_ngram)
eval_dataset.eval_cv(X_no_hla_a2_filter, Y_no_hla_a2_filter)


print 
print "---"
print "No HLA-A2 noisy = positive"
X_no_hla_a2_positive, Y_no_hla_a2_positive = iedb.load_dataset(
                 noisy_labels = 'positive',
                 human = True, 
                 hla_type1 = True,
                 exclude_hla_a2 = True, 
                 only_hla_a2 = False,
                 max_ngram = max_ngram)
eval_dataset.eval_cv(X_no_hla_a2_positive, Y_no_hla_a2_positive)



print 
print "---"
print "No HLA-A2 noisy = negative"
X_no_hla_a2_negtive, Y_no_hla_a2_negative = iedb.load_dataset(
                 noisy_labels = 'negative',
                 human = True, 
                 hla_type1 = True,
                 exclude_hla_a2 = True, 
                 only_hla_a2 = False,
                 max_ngram = max_ngram)
eval_dataset.eval_cv(X_no_hla_a2_positive, Y_no_hla_a2_positive)


print 
print "---"
print "Cross-accuracy for HLA-A2 data"
X_hla_a2, Y_hla_a2 = iedb.load_dataset(
                 noisy_labels = 'keep',
                 human = True, 
                 hla_type1 = True,
                 exclude_hla_a2 = False, 
                 only_hla_a2 = True,
                 max_ngram = max_ngram)
eval_dataset.eval_split(X_no_hla_a2, Y_no_hla_a2, X_hla_a2, Y_hla_a2)


print 
print "---"
print "Cross-accuracy for HLA-A2 data filtered"
X_hla_a2_filtered, Y_hla_a2_filtered = iedb.load_dataset(
                 noisy_labels = 'drop',
                 human = True, 
                 hla_type1 = True,
                 exclude_hla_a2 = False, 
                 only_hla_a2 = True,
                 max_ngram = max_ngram)
eval_dataset.eval_split(X_no_hla_a2_filter, Y_no_hla_a2_filter, X_hla_a2_filtered, Y_hla_a2_filtered)



print 
print "---"
print "Cross-accuracy for HLA-A2 data noisy = positive"
X_hla_a2_positive, Y_hla_a2_positive = iedb.load_dataset(
                 noisy_labels = 'positive',
                 human = True, 
                 hla_type1 = True,
                 exclude_hla_a2 = False, 
                 only_hla_a2 = True,
                 max_ngram = max_ngram)
eval_dataset.eval_split(X_no_hla_a2_positive, Y_no_hla_a2_positive, X_hla_a2_positive, Y_hla_a2_positive)



print 
print "---"
print "Cross-accuracy for HLA-A2 data filtered (assay_group = cytotoxity)"
X_no_hla_a2_cytotoxicity, Y_no_hla_a2_cytotoxicity = iedb.load_dataset(
                 noisy_labels = 'drop',
                 assay_group = 'cytotoxicity', 
                 human = True, 
                 hla_type1 = True,
                 exclude_hla_a2 = True, 
                 only_hla_a2 = False,
                 max_ngram = max_ngram)
                 
X_hla_a2_cytotoxicity, Y_hla_a2_cytotoxicity = iedb.load_dataset(
                 noisy_labels = 'drop',
                 assay_group = 'cytotoxicity', 
                 human = True, 
                 hla_type1 = True,
                 exclude_hla_a2 = False, 
                 only_hla_a2 = True,
                 max_ngram = max_ngram)
                 
eval_dataset.eval_split(X_no_hla_a2_cytotoxicity, Y_no_hla_a2_cytotoxicity, X_hla_a2_cytotoxicity, Y_hla_a2_cytotoxicity)


print 
print "---"
print "Cross-accuracy for HLA-A2 data (noisy = positive, assay_group = cytotoxity)"
X_no_hla_a2_cytotoxicity, Y_no_hla_a2_cytotoxicity = iedb.load_dataset(
                 noisy_labels = 'positive',
                 assay_group = 'cytotoxicity', 
                 human = True, 
                 hla_type1 = True,
                 exclude_hla_a2 = True, 
                 only_hla_a2 = False,
                 max_ngram = max_ngram)
                 
X_hla_a2_cytotoxicity, Y_hla_a2_cytotoxicity = iedb.load_dataset(
                 noisy_labels = 'positive',
                 assay_group = 'cytotoxicity', 
                 human = True, 
                 hla_type1 = True,
                 exclude_hla_a2 = False, 
                 only_hla_a2 = True,
                 max_ngram = max_ngram)
                 
eval_dataset.eval_split(X_no_hla_a2_cytotoxicity, Y_no_hla_a2_cytotoxicity, X_hla_a2_cytotoxicity, Y_hla_a2_cytotoxicity)