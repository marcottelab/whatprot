# -*- coding: utf-8 -*-
"""
@author: Matthew B. Smith (UT Austin)
"""

import copy
import numpy as np
import random

class DyeSeq:
    def __init__(self, n_channels, str_dye_seq, dye_seq_id = -1,
                 src_peptides = []):
        self.n_channels = n_channels
        self.dye_seq = list(str_dye_seq)
        self.dye_seq.reverse()
        self.dye_seq_id = dye_seq_id
        self.src_peptides = src_peptides

    def __getitem__(self, key):
        # Indexing is reversed AND STARTS AT 1. This is not a mistake.
        if key > len(self.dye_seq):
            return '.'
        else:
            return self.dye_seq[-key]
        
    def __setitem__(self, key, value):
        self.dye_seq[-key] = value
        
    def length(self):
        return len(self.dye_seq)
    
    def copy(self):
        return copy.deepcopy(self)
    
    def class_index(self):
        return self.dye_seq_id
    
    def remove_with_prob(self, prob):
        for i in range(len(self.dye_seq)):
            if self.dye_seq[i] == '.':
                continue
            if random.uniform(0, 1) < prob:
                self.dye_seq[i] = '.'
                
    def edman_cycle_with_eff(self, edman_eff):
        if len(self.dye_seq) > 0:
            if random.uniform(0, 1) < edman_eff:
                self.dye_seq.pop(-1)
                
    def get_channel_counts(self):
        counts = [0] * self.n_channels
        for i in range(len(self.dye_seq)):
            if self.dye_seq[i] != '.':
                counts[int(self.dye_seq[i])] += 1
        return counts
    
    def get_channel_intensities(self, mu, sigma, bg_lambda):
        intensities = [0.0] * self.n_channels
        for i in range(len(self.dye_seq)):
            if self.dye_seq[i] != '.':
                intensities[int(self.dye_seq[i])] += np.random.normal(mu,
                            sigma)
        for i in range(len(intensities)):
            if intensities[i] == 0.0 and bg_lambda > 0:
                intensities[i] = np.random.exponential(1 / bg_lambda)
        return intensities
    