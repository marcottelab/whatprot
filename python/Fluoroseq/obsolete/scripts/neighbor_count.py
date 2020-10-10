# -*- coding: utf-8 -*-
"""
@author: Matthew Beauregard Smith (UT Austin)
"""

# This file is meant to be used with a MODIFIED version of the nearest
# neighbors code. It does not work on any .csv file.

from statistics import mean

N_COUNT_FILE = "C:/Users/Matthew/ICES/MarcotteLab/data/classification/control_15_proteins/ann_neighbor_count.csv"
CUTOFF = 100

f = open(N_COUNT_FILE, 'r')
csv = f.readlines()
csv = csv[1:]
neighbor_counts = [0] * len(csv)
for i in range(len(csv)):
    neighbor_counts[i] = int(csv[i].split(",")[1])
f.close()

neighbor_counts = [CUTOFF + (x < CUTOFF) * (x - CUTOFF) for x in neighbor_counts]
print("average number of neighbors: " + str(mean(neighbor_counts)))
