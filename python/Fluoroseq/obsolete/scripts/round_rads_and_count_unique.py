# -*- coding: utf-8 -*-
"""
@author: Matthew Beauregard Smith (UT Austin)
"""

from numpy import around
from numpy import load
from numpy import transpose

NUM_CHANNELS = 3
NUM_CYCLES = 16
MU = 7500.0
RADMAT_FILE = 'C:/Users/Matthew/ICES/MarcotteLab/data/classification/control_15_proteins_full/radmat.npy'

radmat = load(RADMAT_FILE)

radmat = radmat.reshape((radmat.shape[0], NUM_CHANNELS, NUM_CYCLES))
radmat = transpose(radmat, (0, 2, 1))
radmat = radmat.reshape((radmat.shape[0], NUM_CHANNELS * NUM_CYCLES))
radmat = radmat / MU

print("number of radiometries in radmat file: " + str(radmat.shape[0]))

unique_set = set()
for i in range(radmat.shape[0]):
    rad = radmat[i]
    round_rad = around(rad);
    int_rad = (round_rad + .1).astype(int)
    unique_set.add(str(int_rad))

print("number of unique radiometries after rounding to nearest int: " + str(len(unique_set)))
