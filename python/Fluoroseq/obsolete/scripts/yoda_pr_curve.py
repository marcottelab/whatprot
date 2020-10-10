# -*- coding: utf-8 -*-
"""
@author: Matthew Beauregard Smith (UT Austin)
"""

from plotting.plot_pr_curve import plot_pr_curve
from numpy import load

TRUE_Y_FILE = 'C:/Users/Matthew/ICES/MarcotteLab/data/classification/control_15_proteins/true_pep_i.npy'
PREDICTIONS_FILE = 'C:/Users/Matthew/ICES/MarcotteLab/data/classification/control_15_proteins/best_hmm_results.csv'

print("from " + PREDICTIONS_FILE.split("/")[-2] + "/"
              + PREDICTIONS_FILE.split("/")[-1])

true_y = load(TRUE_Y_FILE)

f = open(PREDICTIONS_FILE, "r")
csv = f.readlines()
csv = csv[1:]
predictions = [0] * len(csv)
for i in range(0, len(csv)):
    cells = csv[i].split(",")
    predictions[i] = (int(cells[1]), float(cells[2]))
f.close()

plot_pr_curve(predictions, true_y)
