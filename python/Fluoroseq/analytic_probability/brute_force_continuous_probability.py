# -*- coding: utf-8 -*-
"""
@author: Matthew Beauregard Smith (UT Austin)
"""

from scipy.stats import expon
from scipy.stats import norm

def brute_force_continuous_probability(error_model, dye_seq, observation):
    return brute_force_continuous_probability_starter(error_model, dye_seq,
                                                      observation.intensities)

def brute_force_continuous_probability_starter(error_model, dye_seq,
                                               intensities):
    current_set = [[dye_seq, 1]]
    new_set = []
    for i in range(1, 1 + dye_seq.length()):
        if dye_seq[i] != '.':
            new_set = [0] * 2 * len(current_set)
            for j in range(len(current_set)):
                new_set[2 * j] = current_set[j].copy()
                new_set[2 * j][1] *= (1 - error_model.dud_rate)
                new_set[2 * j + 1] = current_set[j].copy()
                new_set[2 * j + 1][0] = new_set[2 * j + 1][0].copy()
                new_set[2 * j + 1][0][i] = '.'
                new_set[2 * j + 1][1] *= error_model.dud_rate
            current_set = new_set
            new_set = []
    total = 0
    for (ds, dsp) in current_set:
        total += dsp * brute_force_continuous_probability_helper(
                error_model, ds, intensities)
    return total

def brute_force_continuous_probability_helper(error_model, dye_seq,
                                              intensities):
    p = 1
    n_dye = dye_seq.get_channel_counts()
    for c in range(len(n_dye)):
        if n_dye[c] > 0:
            p *= norm.pdf(intensities[0, c], error_model.mu * n_dye[c],
                          error_model.sigma * (n_dye[c] ** .5))
        else:
            p *= expon.pdf(intensities[0, c], 0, 1 / error_model.bg_lambda)
    if len(intensities) == 1:
        return p
    current_set = [[dye_seq, p]]
    new_set = []
    for i in range(1, 1 + dye_seq.length()):
        if dye_seq[i] != '.':
            new_set = [0] * 2 * len(current_set)
            for j in range(len(current_set)):
                new_set[2 * j] = current_set[j].copy()
                new_set[2 * j][1] *= (1 - error_model.bleach_rate)
                new_set[2 * j + 1] = current_set[j].copy()
                new_set[2 * j + 1][0] = new_set[2 * j + 1][0].copy()
                new_set[2 * j + 1][0][i] = '.'
                new_set[2 * j + 1][1] *= error_model.bleach_rate
            current_set = new_set
            new_set = []
    new_set = [0] * 2 * len(current_set)
    for i in range(len(current_set)):
        new_set[2 * i] = current_set[i].copy()
        new_set[2 * i][1] *= (1 - error_model.edman_eff)
        new_set[2 * i + 1] = current_set[i].copy()
        new_set[2 * i + 1][0] = new_set[2 * i + 1][0].copy()
        new_set[2 * i + 1][0].edman_cycle_with_eff(1)
        new_set[2 * i + 1][1] *= error_model.edman_eff
    current_set = new_set
    new_set = []
    total = 0
    for (ds, dsp) in current_set:
        total += dsp * brute_force_continuous_probability_helper(
                error_model, ds, intensities[1:])
    return total
