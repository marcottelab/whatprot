# -*- coding: utf-8 -*-
"""
@author: Matthew Beauregard Smith (UT Austin)
"""

import csv

META_ALL='C:/Users/mbsmi/OneDrive/OdenInstitute/MarcotteLab/data/classification/real-data-4peps-v-50/meta_all.tsv'
RESULTS='C:/Users/mbsmi/OneDrive/OdenInstitute/MarcotteLab/data/classification/real-data-4peps-v-50/rf-predictions.csv'

f_meta = open(META_ALL, 'r')
meta_all = csv.reader(f_meta, delimiter='\t')

f_results = open(RESULTS, 'w')
f_results.write('radmat_iz,best_pep_iz,best_pep_score\n')

idx = 0
for meta_row in meta_all:
    f_results.write(str(idx))
    f_results.write(',')
    f_results.write(meta_row[6])
    f_results.write(',')
    f_results.write(meta_row[2])
    f_results.write('\n')
    idx += 1    

f_meta.close()
f_results.close()
