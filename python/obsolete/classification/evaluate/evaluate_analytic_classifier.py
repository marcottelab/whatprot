# -*- coding: utf-8 -*-
"""
@author: Matthew Beauregard Smith (UT Austin)
"""

from classification.evaluate.plot_pr_curve import plot_pr_curve
from random import sample
from time import process_time

def evaluate_analytic_classifier(classifier, error_model, n_channels,
                                 n_timesteps, dye_seq_list, data, n_dye_seqs,
                                 n_test_per_seq, n_test_limit):
    dye_seqs = [0] * n_dye_seqs
    test = [0] * (n_test_per_seq * n_dye_seqs)
    indices = sample(range(len(data)), n_dye_seqs)
    for i in range(n_dye_seqs):
        dye_seqs[i] = dye_seq_list[indices[i]]
    for i in range(n_dye_seqs):
        test[i * n_test_per_seq : (i + 1) * n_test_per_seq] = (
                sample(data[indices[i]], n_test_per_seq))
        if n_test_limit < n_test_per_seq * n_dye_seqs:
            test = sample(test, n_test_limit)
    print('beginning preparation')
    begin_time = process_time()
    classifier.prepare(error_model, dye_seqs)
    end_time = process_time()
    print('finished preparation')
    print('elapsed time (seconds): ' + str(end_time - begin_time))
    print('time per DyeSeq (seconds): ' + str((end_time - begin_time) /
                            n_dye_seqs))
    print('beginning testing')
    begin_time = process_time()
    predictions = classifier.test(test, True)
    end_time = process_time()
    print('finished testing')
    print('elapsed time (seconds): ' + str(end_time - begin_time))
    print('time per test point (seconds): ' + str((end_time - begin_time) /
                                len(test)))

    plot_pr_curve(predictions, test)