# -*- coding: utf-8 -*-
"""
@author: Matthew Beauregard Smith (UT Austin)
"""

import numpy as np
import random

def bootstrap_radiometries(radmat_file, num_reads_out=-1):
    f = open(radmat_file, 'r')
    line = f.readline()
    num_cycles = int(line)
    line = f.readline()
    num_channels = int(line)
    line = f.readline()
    num_reads = int(line)
    if num_reads_out == -1:
        num_reads_out = num_reads
    f.close()
    radmat = np.genfromtxt(radmat_file, delimiter='\t', dtype=float, skip_header=3)
    for i in range(10):
        output_file = radmat_file.split('.')[0] + '_bootstrap_' + str(i) + '.tsv'
        f = open(output_file, 'w')
        f.write(str(num_cycles) + "\n")
        f.write(str(num_channels) + "\n")
        f.write(str(num_reads_out) + "\n")
        for j in range(num_reads_out):
            x = random.randint(0, num_reads - 1)
            for k in range(num_cycles * num_channels):
                if k > 0:
                    f.write('\t')
                f.write(str(radmat[x,k]))
            f.write('\n')
        f.close()

