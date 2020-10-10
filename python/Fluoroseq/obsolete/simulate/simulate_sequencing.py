# -*- coding: utf-8 -*-
"""
@author: Matthew B. Smith (UT Austin)
"""

from common.channel_intensities import ChannelIntensities
from common.error_model import ErrorModel
from numpy import zeros

def simulate_sequencing(error_model, n_channels, n_timesteps, dye_seq):
    intensities = zeros((n_timesteps, n_channels))
    active_dyes = dye_seq.copy()
    for i in range(n_timesteps):
        if i == 0:
            active_dyes.remove_with_prob(error_model.dud_rate)
        else:
            active_dyes.edman_cycle_with_eff(error_model.edman_eff)
            active_dyes.remove_with_prob(error_model.bleach_rate)
        intensities[i] = active_dyes.get_channel_intensities(error_model.mu,
                   error_model.sigma, error_model.bg_lambda)
    return ChannelIntensities(intensities, dye_seq)
        
def simulate_sequencing_error_free(n_channels, n_timesteps, dye_seq):
    error_model = ErrorModel(1, 0, 0, 1, 0, 0)
    return simulate_sequencing(error_model, n_channels, n_timesteps, dye_seq)