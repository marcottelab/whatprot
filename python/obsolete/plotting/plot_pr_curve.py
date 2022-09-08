# -*- coding: utf-8 -*-
"""
@author: Matthew Beauregard Smith (UT Austin)
"""

import matplotlib.pyplot as plt

def plot_pr_curve(predictions, ground_truth):
    is_correct_and_score = [(False, 0.0)] * len(predictions)

    correct = 0
    incorrect = 0
    for i in range(len(ground_truth)):
        if predictions[i][0] == ground_truth[i]:
            is_correct_and_score[i] = (True, predictions[i][1])
            correct += 1
        else:
            is_correct_and_score[i] = (False, predictions[i][1])
            incorrect += 1
    print('precision at highest recall: ' + str(correct) + '/'
          + str(correct + incorrect) + ', or '
          + str(correct / (correct + incorrect)))

    correct = 0
    incorrect = 0
    total = len(predictions)
    precision = []
    recall = []
    is_correct_and_score.sort(key = lambda entry: entry[1], reverse = True)
    for i in range(len(is_correct_and_score)):
        entry = is_correct_and_score[i]
        next_entry = 0
        if i + 1 < len(is_correct_and_score):
            next_entry = is_correct_and_score[i + 1]
        else:
            next_entry = (False, -1)
        if entry[0] == True:
            correct += 1
        else:
            incorrect += 1
        if entry[1] != next_entry[1]:
            precision.append(correct / (correct + incorrect))
            recall.append(correct / total)
    plt.plot(recall, precision, '-')
    plt.xlabel('recall')
    plt.ylabel('precision')
    plt.xlim(0, 1)
    plt.ylim(0, 1)
    plt.show()
