# -*- coding: utf-8 -*-
"""
@author: Matthew Beauregard Smith (UT Austin)
"""

from common.peptide import Peptide
from simulate.label_peptides import label_peptides

def dye_seqs_from_peptides(peptide_file, label_set, dye_seqs_file):
    f = open(peptide_file, 'r')
    # f.readline()  # header
    # f.readline()  # Zack's null line
    line = "placeholder"
    peptides = []
    pep_id = 0
    while line != '\n' and line != '':
        line = f.readline()[0 : -1]
        peptides += [Peptide(line, pep_id=pep_id)]
        pep_id += 1
    f.close()
    dye_seqs = label_peptides(peptides, label_set)
    
    f = open(dye_seqs_file, 'w')
    f.write(str(len(label_set)) + "\n")  # num channels
    f.write(str(len(dye_seqs)) + "\n")
    for dye_seq in dye_seqs:
        dye_seq.dye_seq.reverse()
        f.write("".join(dye_seq.dye_seq) + "\t")
        f.write(str(len(dye_seq.src_peptides)) + "\t")
        f.write(str(dye_seq.src_peptides[0].pep_id) + "\n")
    f.close()

