# -*- coding: utf-8 -*-
"""
@author: Matthew B. Smith (UT Austin)
"""

from random import sample
from simulate.label_peptides import label_peptides
from simulate.load_proteins import load_proteins
from simulate.simulate_protease import simulate_protease

def generate_dye_seqs(n_dye_seqs, labels_by_channel):
    proteins = load_proteins()
    peptides = simulate_protease(proteins)
    dye_seqs = label_peptides(peptides, labels_by_channel)
    print("total number of possible dye_seqs: " + str(len(dye_seqs)))
    return sample(dye_seqs, n_dye_seqs)