# -*- coding: utf-8 -*-
"""
@author: Matthew Beauregard Smith (UT Austin)
"""

import matplotlib.pyplot as plt
from numpy import load

def plot_pr_curve(predictions_file,
                  true_y_file,
                  directory = "",
                  title = "Precision / Recall"):
    plot_pr_curves([predictions_file], true_y_file, directory, title)

def plot_pr_curves(predictions_files,
                   true_y_file,
                   directory = "",
                   title = "Precision / Recall",
                   labels = []):
    fig, ax = plt.subplots()
    cmloc = [x / len(predictions_files) for x in range(len(predictions_files))]
    ax.set_prop_cycle("color", plt.get_cmap("hsv")(cmloc))
    fig.set_size_inches(8, 8)
    fig.set_dpi(100)
    for i in range(len(predictions_files)):
        predictions_file = predictions_files[i]
        label = None
        if labels != []:
            label = labels[i]
        plot_pr_curve_noshow(ax,
                             predictions_file,
                             true_y_file,
                             directory = directory,
                             label = label)
    ax.set_xlim(left = 0)
    ax.set_ylim(bottom = 0)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["bottom"].set_linewidth(1.5)
    ax.spines["left"].set_linewidth(1.5)
    ax.xaxis.set_tick_params(width=1.5)
    ax.yaxis.set_tick_params(width=1.5)
    ax.set_title(title,
                 fontsize = 25,
                 pad = 15)
    ax.set_xlabel('recall',
                  fontsize = 15,
                  labelpad = 15)
    ax.set_ylabel('precision',
                  fontsize = 15,
                  labelpad = 15)
    if labels != []:
        plt.legend(fontsize = 15)

def plot_pr_curve_noshow(ax,
                         predictions_file,
                         true_y_file,
                         directory = "",
                         label = None):
    true_y = None
    if (true_y_file.split(".")[-1] == "npy"):
        true_y = load(directory + true_y_file)
    else:  # then it's a tsv file
        f = open(directory + true_y_file, "r")
        tsv = f.readlines()
        tsv = tsv[1:]
        true_y = [int(x) for x in tsv]
    f = open(directory + predictions_file, "r")
    csv = f.readlines()
    csv = csv[1:]
    predictions = [0] * len(csv)
    for i in range(0, len(csv)):
        cells = csv[i].split(",")
        predictions[i] = (int(cells[1]), float(cells[2]))
    f.close()
    precision, recall = compute_pr_curve(predictions, true_y)
    ax.plot(recall, precision, '-', label = label, linewidth = 2)

def compute_pr_curve(predictions, ground_truth):
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
#    print('precision at highest recall: ' + str(correct) + '/'
#          + str(correct + incorrect) + ', or '
#          + str(correct / (correct + incorrect)))

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
            if len(precision) == 0:
                precision.append([correct / (correct + incorrect)])
                recall.append(0)
#            if correct / total == recall[-1]:
#                precision[-1].append(correct / (correct + incorrect))
#            else:
            precision.append([correct / (correct + incorrect)])
            recall.append(correct / total)
    for i in range(len(recall)):
        precision[i] = sum(precision[i]) / len(precision[i])
    return (precision, recall)
