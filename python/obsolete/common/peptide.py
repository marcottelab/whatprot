# -*- coding: utf-8 -*-
"""
@author: Matthew B. Smith (UT Austin)
"""

class Peptide:
    def __init__(self, seq, src_proteins=[], pep_id=-1):
        self.seq = seq
        self.src_proteins = src_proteins
        self.pep_id = pep_id
