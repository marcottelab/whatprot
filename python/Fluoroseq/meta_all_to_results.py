# -*- coding: utf-8 -*-
"""
@author: Matthew Beauregard Smith (UT Austin)
"""

import csv

# Converts classification predictions from Erisyon into the .csv format given
# by classification with whatprot.
def meta_all_to_results(meta_all, results):
    f_meta = open(meta_all, 'r')
    meta_all = csv.reader(f_meta, delimiter='\t')

    f_results = open(results, 'w')
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
