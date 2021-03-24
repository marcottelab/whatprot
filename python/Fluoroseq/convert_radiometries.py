# -*- coding: utf-8 -*-
"""
@author: Matthew Beauregard Smith (UT Austin)
"""

from numpy import load
from numpy import transpose

NUM_CHANNELS = 2
NUM_CYCLES = 16
BETA = 15000
RADMAT_FILE = 'C:/Users/Matthew/OneDrive/OdenInstitute/MarcotteLab/data/classification/yoda_2_decoys_lysc_de_c_15k/radmat.npy'
OUTPUT_FILE = 'C:/Users/Matthew/OneDrive/OdenInstitute/MarcotteLab/data/classification/yoda_2_decoys_lysc_de_c_15k/radiometries.tsv'

radmat = load(RADMAT_FILE)

radmat = transpose(radmat, (0, 2, 1))

# Fix intensities
radmat = radmat / BETA

f = open(OUTPUT_FILE, 'w')
f.write(str(NUM_CYCLES) + "\n")
f.write(str(NUM_CHANNELS) + "\n")
f.write(str(radmat.shape[0]) + "\n")
for rad in radmat:
    for i in range(NUM_CYCLES):
        for j in range(NUM_CHANNELS):
            if i > 0 or j > 0:
                f.write("\t")
            # Zero out anything less than half the one dye intensity
            if rad[i, j] < 0.5:
                f.write(str(0.0))
            else:
                f.write(str(rad[i, j]))
    f.write("\n")
f.close()
