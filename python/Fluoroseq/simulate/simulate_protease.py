# -*- coding: utf-8 -*-
"""
@author: Matthew B. Smith (UT Austin)
"""

# Currently using protease lysc only. Clips after K (Lysine).

from collections import defaultdict
from common.peptide import Peptide

def simulate_protease(proteins):
    peptide_dict = defaultdict(lambda: [])
    for protein in proteins:
        peptides = protein.seq.split("K")
        for i in range(len(peptides) - 1):
            peptides[i] += "K"
        if peptides[-1] == "":
            peptides = peptides[0:-1]
        for peptide in peptides:
            peptide_dict[peptide] += [protein.name]
    peptides = [0] * len(peptide_dict)
    index = 0
    for peptide in peptide_dict:
        peptides[index] = Peptide(peptide, peptide_dict[peptide])
        index += 1
    return peptides