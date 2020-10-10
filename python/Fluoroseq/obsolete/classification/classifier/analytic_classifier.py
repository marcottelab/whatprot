# -*- coding: utf-8 -*-
"""
@author: Matthew Beauregard Smith (UT Austin)
"""

from analytic_probability.continuous_probability import continuous_probability

class AnalyticClassifier:
    def __init__(self):
        self.error_model = 0
        self.dye_seq_list = 0
    
    def prepare(self, error_model, dye_seq_list):
        self.error_model = error_model
        self.dye_seq_list = dye_seq_list
        
    def test(self, test_data, include_scores = False):
        results = [0] * len(test_data)
        for i in range(len(test_data)):
            results[i] = self.classify(test_data[i], include_scores)
        return results
    
    def classify(self, test, include_scores = False):
        best_score = 0
        best_dye_seq = -1
        total_score = 0
        for dye_seq in self.dye_seq_list:
            score = continuous_probability(self.error_model, dye_seq, test)
            total_score += score
            if score > best_score:
                best_score = score
                best_dye_seq = dye_seq.dye_seq_id
        if include_scores:
            result = (best_dye_seq, best_score / total_score)
        else:
            result = best_dye_seq
        return result
    