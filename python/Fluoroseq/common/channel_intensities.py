# -*- coding: utf-8 -*-
"""
@author: Matthew B. Smith (UT Austin)
"""

class ChannelIntensities:
    def __init__(self, intensities, src_dye_seq=None):
        self.intensities = intensities
        self.src_dye_seq = src_dye_seq
    
    def class_index(self):
        return self.src_dye_seq.dye_seq_id
    
    def feature_vector(self):
        return self.intensities.flatten()
    
    def __str__(self):
        return self.intensities.str()
