# -*- coding: utf-8 -*-
"""
@author: Matthew B. Smith (UT Austin)
"""

from simulate.simulate_sequencing import simulate_sequencing

def generate_sequences(error_model, n_channels, n_timesteps, dye_seqs,
                       n_samples_per_dye_seq):
    result = [0] * len(dye_seqs)
    for i in range(len(dye_seqs)):
        result[i] = [0] * n_samples_per_dye_seq
        for j in range(n_samples_per_dye_seq):
            result[i][j] = simulate_sequencing(error_model, n_channels,
                  n_timesteps, dye_seqs[i])
    return result