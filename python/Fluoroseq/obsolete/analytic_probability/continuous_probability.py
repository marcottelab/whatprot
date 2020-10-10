# -*- coding: utf-8 -*-
"""
@author: Matthew Beauregard Smith (UT Austin)
"""

from functools import lru_cache
from numpy import array
from numpy import log
from numpy import sum as ndsum
from numpy import tensordot
from numpy import zeros
from scipy.sparse import diags
from scipy.stats import binom
from scipy.stats import expon
from scipy.stats import norm
from simulate.simulate_sequencing import simulate_sequencing_error_free

def continuous_probability(error_model, dye_seq, observation):
    n_channels = len(observation.intensities[0])
    n_timesteps = len(observation.intensities)
    # n_dye: number of physical dyes remaining, indexed by [edman cycle,
    # channel]
    n_dye = array(
            simulate_sequencing_error_free(
                    n_channels, n_timesteps, dye_seq).intensities,
            dtype=int)
    intensities = observation.intensities
    # p: probability density, indexed by [edman cycle, (number of remaining
    # active dyes in each channel ...)]
    p = zeros((n_timesteps,) + tuple([n + 1 for n in n_dye[0]]))
    p[0][tuple(n_dye[0])] = 1
    dud_tensor = build_dye_removal_tensor(max(n_dye[0]),
                                          1 - error_model.dud_rate)
    bleach_tensor = build_dye_removal_tensor(max(n_dye[0]),
                                             1 - error_model.bleach_rate)
    # Adjust initial probabilities accounting for dud rates.
    remove_dye(p, 1, n_channels, dud_tensor)
    # Adjust initial probabilities for initial observation.
    if error_model.bg_lambda > 0:
        observe_intensities(p[0:1], intensities[0], error_model)
    else:
        safe_observe_intensities(p[0:1], intensities[0], error_model)
    # t: timestep
    for t in range(1, n_timesteps):
        # Account for photobleaching.
        remove_dye(p, t, n_channels, bleach_tensor)
        # Account for detach rate.
        detach_peptide(p[0:t], error_model.detach_rate)
        # Account for edman cycle.
        edman_cycle(p[0:(t+1)], dye_seq, n_dye, error_model.edman_eff)
        # Adjust for observation.
        if error_model.bg_lambda > 0:
            observe_intensities(p[0:(t+1)], intensities[t], error_model)
        else:
            safe_observe_intensities(p[0:(t+1)], intensities[t], error_model)
    # Probability is sum of probability over all amounts of edman-cycles and
    # dye-counts.
    return ndsum(p)

def build_dye_removal_tensor(size, p):
    tensor = zeros((size + 1, size + 1))
    for k in range(size + 1):
        for n in range(k, size + 1):
            tensor[k, n] = fast_binom(k, n, p)
    return tensor

def remove_dye(p, t, n_channels, dud_tensor):
    for c in range(n_channels):
        lenc = p.shape[1 + c]
        s = list(range(len(p.shape)))
        s[1 + c] = 0
        s[0] = 1 + c
        q = p[0:t]
        q = q.transpose(s)
        q = tensordot(dud_tensor[0:lenc, 0:lenc], q, ([1], [0]))
        q = q.transpose(s)
        p[0:t] = q

def detach_peptide(p, detach_rate):
    total = ndsum(p)
    p *= (1 - detach_rate)
    p[tuple(zeros(len(p.shape), int))] += detach_rate * total


# this function currently doesn't support lognormal. Must use
# safe_observe_intensities instead.
def observe_intensities(p, intensities, error_model):
    for c in range(len(p.shape) - 1):
        lenc = p.shape[1 + c]
        normidx = array(range(1, lenc))
        pdf_tensor = zeros(lenc)
        pdf_tensor[0] = expon.pdf(intensities[c], 0, 1 / error_model.bg_lambda)
        pdf_tensor[1:lenc] = norm.pdf(intensities[c],
                                      error_model.mu * normidx,
                                      error_model.sigma * (normidx ** .5))
        multiply(pdf_tensor, p, 1 + c)

def safe_observe_intensities(p, intensities, error_model):
    for c in range(len(p.shape) - 1):
        lenc = p.shape[1 + c]
        normidx = array(range(1, lenc))
        pdf_tensor = zeros(lenc)
        if intensities[c] == 0:
            pdf_tensor[0] = 1
            pdf_tensor[1:lenc] = 0
        else:
            pdf_tensor[0] = 0
            if error_model.is_lognormal:
                pdf_tensor[1:lenc] = my_lognorm(
                                          x = intensities[c],
                                          mu = error_model.log_mu
                                               + log(normidx),
                                          sigma = error_model.sigma)
            else:
                pdf_tensor[1:lenc] = norm.pdf(
                                          x = intensities[c],
                                          loc = error_model.mu * normidx,
                                          scale = error_model.sigma
                                                  * (normidx ** .5))
        multiply(pdf_tensor, p, 1 + c)

def edman_cycle(p, dye_seq, n_dye, edman_eff):
    p_edman = p[:-1] * edman_eff
    for e in range(len(p_edman)):
        if dye_seq[e + 1] == '.':
            continue
        c = int(dye_seq[e + 1])
        lenc = p.shape[1 + c]
        n = n_dye[e, c]
        d1 = zeros(lenc)
        d1[:(n+1)] = [i / n for i in reversed(range(n + 1))]
        d2 = zeros(lenc - 1)
        d2[:(n)] = [i / n for i in range(1, n + 1)]
        d = diags([d1, d2], [0, 1]).toarray()
        p_edman[e] = matrix_multiply(d, p_edman[e], c)
    p[:] *= (1 - edman_eff)
    p[1:] += p_edman

def multiply(a, b, axis):
    t = list(range(len(b.shape)))
    t[axis] = t[-1]
    t[-1] = axis
    b = b.transpose(t)
    b *= a
    b = b.transpose(t)

def matrix_multiply(a, b, axis):
    t = list(range(len(b.shape)))
    t[axis] = 0
    t[0] = axis
    b = b.transpose(t)
    b = tensordot(a, b, ([1], [0]))
    b = b.transpose(t)
    return b

def my_lognorm(x, mu, sigma):
    return norm.pdf(log(x), mu, sigma) / x

@lru_cache(maxsize = None)
def fast_binom(k, n, p):
    return binom.pmf(k, n, p)
