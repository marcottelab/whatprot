# -*- coding: utf-8 -*-
"""
@author: Matthew B. Smith (UT Austin)
"""

from math import log

class ErrorModel:
    def __init__(self,
                 edman_eff,
                 bleach_rate,
                 dud_rate,
                 mu,
                 sigma,
                 bg_lambda,
                 detach_rate = 0.0,
                 is_lognormal = False):
        self.edman_eff = edman_eff
        self.bleach_rate = bleach_rate
        self.dud_rate = dud_rate
        self.mu = mu
        self.log_mu = log(mu)
        self.sigma = sigma
        self.bg_lambda = bg_lambda
        self.detach_rate = detach_rate
        self.is_lognormal = is_lognormal