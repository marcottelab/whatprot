# -*- coding: utf-8 -*-
"""
@author: Matthew Beauregard Smith (UT Austin)
"""

from common.peptide import Peptide
from plotting.plot_pr_curve import plot_pr_curve
from numpy import load
from simulate.label_peptides import label_peptides


TRUE_Y_FILE = 'C:/Users/Matthew/ICES/MarcotteLab/data/classification/control_15_proteins/true_pep_i.npy'
NUM_PEPTIDES = 705
NUM_CHANNELS = 3
LABEL_SET = ['DE','Y','C']
PEPTIDE_FILE = 'C:/Users/Matthew/ICES/MarcotteLab/data/classification/control_15_proteins/peps.csv'

true_y = load(TRUE_Y_FILE)
class GroundTruth:
    def __init__(self, value):
        self.value = value
    def class_index(self):
        return self.value
ground_truth = [0] * len(true_y)
for i in range(0, len(true_y)):
    ground_truth[i] = GroundTruth(true_y[i])

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
id_to_prediction = {}
for dye_seq in dye_seqs:
    for peptide in dye_seq.src_peptides:
        id_to_prediction[int(peptide.pep_id)] = (
                int(dye_seq.src_peptides[0].pep_id),
                1 / len(dye_seq.src_peptides))
predictions = [0] * len(ground_truth)
for i in range(len(ground_truth)):
    predictions[i] = id_to_prediction[ground_truth[i].value]

plot_pr_curve(predictions, ground_truth)
