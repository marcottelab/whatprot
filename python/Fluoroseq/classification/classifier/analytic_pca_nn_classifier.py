# -*- coding: utf-8 -*-
"""
@author: Matthew Beauregard Smith (UT Austin)
"""

from analytic_probability.continuous_probability import continuous_probability
from collections import defaultdict
from heapq import nlargest
from scipy.spatial import KDTree
from sklearn.decomposition import PCA

class AnalyticPCANNClassifier:
    def __init__(self, nn_k, d, leafsize, analytic_k):
        self.nn_k = nn_k
        self.d = d
        self.pca = PCA(n_components = d)
        self.leafsize = leafsize
        self.analytic_k = analytic_k
        self.training_data = 0
        self.reduced_training_data = 0
        self.nn_kdtree = 0
        self.error_model = 0
        self.dye_seq_dict = 0
    
    def prepare(self, error_model, dye_seq_list):
        self.error_model = error_model
        self.dye_seq_dict = {}
        for dye_seq in dye_seq_list:
            self.dye_seq_dict[dye_seq.class_index()] = dye_seq
        
    def train(self, training_data):
        self.training_data = training_data
        feature_vectors = [0] * len(training_data)
        for i in range(len(training_data)):
            feature_vectors[i] = training_data[i].feature_vector()
        self.reduced_training_data = self.pca.fit_transform(feature_vectors)
        self.nn_kdtree = KDTree(self.reduced_training_data, self.leafsize)
    
    def test(self, test_data, include_scores = False):
        results = [0] * len(test_data)
        for i in range(len(test_data)):
            results[i] = self.classify(test_data[i], include_scores)
        return results
    
    def classify(self, test, include_scores = False):
        reduced_test = self.pca.transform(test.feature_vector().reshape(1, -1))
        reduced_test = reduced_test.flatten()
        _, indices = self.nn_kdtree.query(reduced_test, self.nn_k)
        if self.nn_k == 1:
            return self.training_data[indices]
        nn_scores = defaultdict(lambda: 0)
        for i in range(self.nn_k):
            train = self.training_data[indices[i]]
            nn_scores[train.class_index()] += 1
        nn_k_best = []
        if self.analytic_k < len(nn_scores):
            nn_k_best = nlargest(self.analytic_k, nn_scores,
                                 key = nn_scores.get)
        else:
            nn_k_best = nn_scores.keys()
        best_analytic_score = 0
        best_analytic_dye_seq = -1
        total_analytic_score = 0
        for i in nn_k_best:
            dye_seq = self.dye_seq_dict[i]
            analytic_score = continuous_probability(self.error_model, dye_seq,
                                                    test)
            total_analytic_score += analytic_score
            if analytic_score > best_analytic_score:
                best_analytic_score = analytic_score
                best_analytic_dye_seq = dye_seq.dye_seq_id
        if include_scores:
            result = (best_analytic_dye_seq,
                      best_analytic_score / total_analytic_score)
        else:
            result = best_analytic_dye_seq
        return result
    