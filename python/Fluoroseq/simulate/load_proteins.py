# -*- coding: utf-8 -*-
"""
@author: Matthew B. Smith (UT Austin)
"""

from common.protein import Protein

def load_proteins():
    f = open("simulate\\human_proteome\\UP000005640_9606.fasta", "r")
    aa_seqs = {}
    if f.mode == 'r':
        protein_id = ""
        while True:
            line = f.readline()
            if line == '':
                break
            if line[0] == '>':
                protein_id = line.split("|")[1]
                aa_seqs[protein_id] = ""
            else:
                aa_seqs[protein_id] += line[0:-1]
    f.close()
    proteins = [0] * len(aa_seqs)
    index = 0
    for protein_id in aa_seqs:
        proteins[index] = Protein(protein_id, aa_seqs[protein_id])
        index += 1
    return proteins