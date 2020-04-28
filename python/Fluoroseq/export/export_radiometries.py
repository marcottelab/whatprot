# -*- coding: utf-8 -*-
"""
@author: Matthew Beauregard Smith (UT Austin)
"""

from numpy import load
from numpy import transpose

NUM_CHANNELS = 3
NUM_CYCLES = 16
MU = 7500.0
RADMAT_FILE = 'analytic_probability/control_15_proteins/radmat.npy'
OUTPUT_FILE = 'analytic_probability/control_15_proteins/radiometries.tsv'

radmat = load(RADMAT_FILE)

radmat = radmat.reshape((radmat.shape[0], NUM_CHANNELS, NUM_CYCLES))
radmat = transpose(radmat, (0, 2, 1))
radmat = radmat.reshape((radmat.shape[0], NUM_CHANNELS * NUM_CYCLES))
radmat = radmat / MU

f = open(OUTPUT_FILE, 'w')
f.write(str(NUM_CYCLES) + "\n")
f.write(str(NUM_CHANNELS) + "\n")
f.write(str(radmat.shape[0]) + "\n")
for rad in radmat:
    f.write(str(rad[0]))
    for i in range(1, NUM_CHANNELS * NUM_CYCLES):
        f.write("\t" + str(rad[i]))
    f.write("\n")
f.close()
