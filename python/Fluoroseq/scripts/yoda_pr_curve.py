# -*- coding: utf-8 -*-
"""
@author: Matthew Beauregard Smith (UT Austin)
"""

from plotting.plot_pr_curve import plot_pr_curve
from numpy import load

TRUE_Y_FILE = 'C:/Users/Matthew/ICES/MarcotteLab/data/classification/control_3_proteins/true_pep_i.npy'
PREDICTIONS_FILE = 'C:/Users/Matthew/ICES/MarcotteLab/data/classification/control_3_proteins/hybrid_test.csv'

true_y = load(TRUE_Y_FILE)
class GroundTruth:
    def __init__(self, value):
        self.value = value
    def class_index(self):
        return self.value
ground_truth = [0] * len(true_y)
for i in range(0, len(true_y)):
    ground_truth[i] = GroundTruth(true_y[i])

f = open(PREDICTIONS_FILE, "r")
csv = f.readlines()
csv = csv[1:]
predictions = [0] * len(csv)
for i in range(0, len(csv)):
    cells = csv[i].split(",")
    predictions[i] = (int(cells[1]), float(cells[2]))
f.close()

plot_pr_curve(predictions, ground_truth)
