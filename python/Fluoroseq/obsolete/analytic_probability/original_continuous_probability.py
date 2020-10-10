# -*- coding: utf-8 -*-
"""
@author: Matthew Beauregard Smith (UT Austin)
"""

from functools import lru_cache
from numpy import array
from numpy import sum as ndsum
from numpy import zeros
from scipy.stats import binom
from scipy.stats import expon
from scipy.stats import norm
from simulate.simulate_sequencing import simulate_sequencing_error_free

def original_continuous_probability(error_model, dye_seq, observation):
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
    # Iteration to compute initial probabilities from initial intensity and
    # the dud rate.
    # x: Iterates number of active (non-dud) dyes, indexed by channel.
    for x in box_range(n_dye[0]):
        p[0][x] = 1
        # For each channel, compute the probability of our initial observation
        # (and take product for full probability density)
        for c in range(n_channels):
            # Binomial to compute probability of this particular number of dud-
            # dyes.
            p[0][x] *= fast_binom(x[c], n_dye[0, c], 1 - error_model.dud_rate)
            # Model for density of intensity.
            if x[c] > 0:
                p[0][x] *= norm.pdf(intensities[0, c], x[c] * error_model.mu,
                                    x[c] ** 0.5 * error_model.sigma)
            else:
                p[0][x] *= expon.pdf(intensities[0, c], 0,
                                     1 / error_model.bg_lambda)
    # t: Iterates timesteps (should be used when indexing intensities)
    for t in range(1, n_timesteps):
        # e: Iterates edman cycles (should be used when indexing dye counts).
        # Is reversed to keep data from being overwritten as we advance.
        for e in reversed(range(t + 1)):
            # x: Iterates dye counts (indexed by channel) to fill in.
            for x in box_range(n_dye[e]):
                # y: Iterates dye counts (indexed by channel) of the previous
                # timestep. Relation to x is important; iterating in this
                # particular manner guarantees that x[c] <= y[c] for all c,
                # so we can overwrite data in p. This would NOT be true if we
                # swapped the x and y for loops.
                for y in box_range(x, n_dye[max(0, e - 1)]):
                    p_e_succeed = 0  # Probability with edman success.
                    # Conditional: Must have had at least one successful edman
                    # cycle to consider successful edman cycle.
                    if e > 0:
                        # Relevant probability, times chance of edman success.
                        p_e_succeed = p[e - 1][y] * error_model.edman_eff
                        for c in range(n_channels):
                            # Conditional: If dye passed was the same as the
                            # channel we are on, we must consider the two cases
                            # with and without the dye being active (not a dud
                            # or photo-bleached).
                            if str(c) == dye_seq[e]:
                                p_with_dye = 0
                                # Condition: Must have had at least one active
                                # dye last cycle to decrease intensity from
                                # edman cycle.
                                if y[c] > 0:
                                    # First, consider probability dye is
                                    # active.
                                    p_with_dye = (p_e_succeed * y[c]
                                                  / n_dye[e - 1, c])
                                    # Second, consider probability of some
                                    # number of photo-bleaching events
                                    # occurring.
                                    p_with_dye *= fast_binom(
                                            x[c], y[c] - 1,
                                            1 - error_model.bleach_rate)
                                # First, consider probability dye is NOT
                                # active.
                                p_without_dye = (p_e_succeed
                                                 * (n_dye[e - 1, c] - y[c])
                                                 / n_dye[e - 1, c])
                                # Second, consider probability of some number
                                # of photo-bleaching events occurring.
                                p_without_dye *= fast_binom(
                                        x[c], y[c],
                                        1 - error_model.bleach_rate)
                                # Combine two scenarios to get total
                                # probability with successful edman cycle.
                                p_e_succeed = p_with_dye + p_without_dye
                            else:
                                # The else case is much simpler. Just need to
                                # compute probability of photo-bleaching
                                # events.
                                p_e_succeed *= fast_binom(
                                        x[c], y[c],
                                        1 - error_model.bleach_rate)
                    p_e_fail = 0  # probability with edman failure.
                    # Conditional: Must exclude case with same number of edman
                    # cycles and time-steps. Then there can not have been an
                    # edman failure. Additionally, other things in this block
                    # require numbers of dyes from previous edman cycle, but
                    # since we are considering the case where the edman cycle
                    # failed, we need to make sure we are within the number of
                    # dyes of the current edman cycle, NOT the previous.
                    if e < t and y in box_range(x, n_dye[e]):
                        # Relevant probability, times chance of edman failure.
                        p_e_fail = p[e][y] * (1 - error_model.edman_eff)
                        for c in range(n_channels):
                            # Probability of exact number of photo-bleaching
                            # events.
                            p_e_fail *= fast_binom(
                                    x[c], y[c], 1 - error_model.bleach_rate)
                    # y = x is the first iteration of the y loop. We want to
                    # accumulate all the relevant probabilities over y, so in
                    # the first iteration we must set our result to 0.
                    if x == y:
                        p[e][x] = 0
                    # Total probability is sum of probability with and without
                    # edman failure.
                    p[e][x] += p_e_succeed + p_e_fail
                for c in range(n_channels):
                    # Model for density of intensity.
                    if x[c] > 0:
                        p[e][x] *= norm.pdf(intensities[t, c],
                                            x[c] * error_model.mu,
                                            x[c] ** .5 * error_model.sigma)
                    else:
                        p[e][x] *= expon.pdf(intensities[t, c], 0,
                                             1 / error_model.bg_lambda)
    # Probability is sum of probability over all amounts of edman-cycles and
    # dye-counts.
    return ndsum(p)

# Iterates through a multi-dimensional box.
class BoxRange:
    def __init__(self, start, stop):
        self.start = start
        self.stop = stop
        self.current = list(start)
        self.active = False
    
    def __iter__(self):
        return self
    
    def __next__(self):
        if not self.active:
            self.active = True
            return tuple(self.current)
        o = len(self.current) - 1
        while True:
            self.current[o] += 1
            if self.current[o] <= self.stop[o]:
                break
            self.current[o] = self.start[o]
            o -= 1
            if o < 0:
                raise StopIteration
        return tuple(self.current)
    
    def __contains__(self, item):
        for o in range(len(self.stop)):
            if not (self.start[o] <= item[o] <= self.stop[o]):
                return False
        return True

def box_range(shape1, shape2 = ()):
    if shape2 == ():
        return BoxRange(tuple([0] * len(shape1)), shape1)
    else:
        return BoxRange(shape1, shape2)
    
@lru_cache(maxsize = None)
def fast_binom(k, n, p):
    return binom.pmf(k, n, p)