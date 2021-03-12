# -*- coding: utf-8 -*-
"""
@author: Matthew B. Smith (UT Austin)
"""

class DyeTrack:
    def __init__(self, intensities, src_dye_seq=None):
        self.intensities = intensities
        self.src_dye_seq = src_dye_seq

    def feature_vector(self):
        return self.intensities.flatten()

    def __str__(self):
        return str(self.intensities)
