# -*- coding: utf-8 -*-
"""
@author: Matthew Beauregard Smith (UT Austin)
"""

from numpy import load
from numpy import transpose
from numpy import genfromtxt
from numpy import reshape

NUM_CHANNELS = 2
NUM_MOCKS = 0
NUM_CYCLES = 16
BETA = 15000
# RADMAT_FILE = 'C:/Users/mbsmi/OneDrive/OdenInstitute/MarcotteLab/data/classification/jim_jhm2022_07_04a_01_jsp127_1p2m10e_640_startcycle2/full_signal_radmat.npy'
TSV_RADMAT_FILE = 'C:/Users/mbsmi/OneDrive/OdenInstitute/MarcotteLab/data/classification/real-data-4peps-v-50/vector_all.tsv'
OUTPUT_FILE = 'C:/Users/mbsmi/OneDrive/OdenInstitute/MarcotteLab/data/classification/real-data-4peps-v-50/radiometries.tsv'

# radmat = load(RADMAT_FILE)

radmat = genfromtxt(TSV_RADMAT_FILE, delimiter='\t', dtype=float)[:,1:]
radmat = reshape(radmat, (radmat.shape[0], NUM_CHANNELS, NUM_MOCKS + NUM_CYCLES))

radmat = transpose(radmat, (0, 2, 1))

# Fix intensities
radmat = radmat / BETA

f = open(OUTPUT_FILE, 'w')
f.write(str(NUM_CYCLES) + "\n")
f.write(str(NUM_CHANNELS) + "\n")
f.write(str(radmat.shape[0]) + "\n")
for rad in radmat:
    for i in range(NUM_MOCKS, NUM_MOCKS + NUM_CYCLES):
        for j in range(NUM_CHANNELS):
            if i > 0 or j > 0:
                f.write("\t")
            f.write(str(rad[i, j]))
    f.write("\n")
f.close()
