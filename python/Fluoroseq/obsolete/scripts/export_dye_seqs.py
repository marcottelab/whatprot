# -*- coding: utf-8 -*-
"""
@author: Matthew Beauregard Smith (UT Austin)
"""

from common.peptide import Peptide
from simulate.label_peptides import label_peptides

NUM_PEPTIDES = 70263
NUM_CHANNELS = 3
LABEL_SET = ['DE','Y','C']
PEPTIDE_FILE = 'C:/Users/Matthew/ICES/MarcotteLab/data/classification/n2000_proteins_trypsin/human_random_2000_proteins_trypsin.csv'
OUTPUT_FILE = 'C:/Users/Matthew/ICES/MarcotteLab/data/classification/n2000_proteins_trypsin/dye_seqs.tsv'

f = open(PEPTIDE_FILE, 'r')
f.readline()  # header
f.readline()  # Zack's null line
line = f.readline()
peptides = [0] * NUM_PEPTIDES
i = 0
while line != '\n' and line != '':
    items = line.split(",")
    pep_id = items[0]
    pep_str = items[-1]
    peptides[i] = Peptide(pep_str, pep_id=pep_id)
    line = f.readline()
    i += 1
f.close()
dye_seqs = label_peptides(peptides, LABEL_SET)

f = open(OUTPUT_FILE, 'w')
f.write(str(NUM_CHANNELS) + "\n")
f.write(str(len(dye_seqs)) + "\n")
for dye_seq in dye_seqs:
    dye_seq.dye_seq.reverse()
    f.write("".join(dye_seq.dye_seq) + "\t")
    f.write(str(len(dye_seq.src_peptides)) + "\t")
    f.write(str(dye_seq.src_peptides[0].pep_id) + "\n")
f.close()

