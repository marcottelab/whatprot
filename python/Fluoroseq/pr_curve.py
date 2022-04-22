# -*- coding: utf-8 -*-
"""
@author: Matthew Beauregard Smith (UT Austin)
"""

import matplotlib.pyplot as plt
from numpy import load

plt.rcParams["font.family"] = "Calibri"

def plot_pr_curve(predictions_file,
                  true_y_file,
                  dye_seqs_file,
                  directory = "",
                  title = "Precision / Recall"):
    plot_pr_curves([predictions_file],
                   true_y_file,
                   dye_seqs_file,
                   directory,
                   title)

def plot_pr_curves(predictions_files,
                   true_y_file,
                   dye_seqs_file,
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
        cur_true_y_file = ''
        if isinstance(true_y_file, str):
            cur_true_y_file = true_y_file
        else:
            cur_true_y_file = true_y_file[i]
        cur_dye_seqs_file = ''
        if isinstance(dye_seqs_file, str):
            cur_dye_seqs_file = dye_seqs_file
        else:
            cur_dye_seqs_file = dye_seqs_file[i]
        plot_pr_curve_noshow(ax,
                                predictions_file,
                                cur_true_y_file,
                                cur_dye_seqs_file,
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
                         dye_seqs_file,
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
        f.close()
    f = open(directory + dye_seqs_file, "r")
    f.readline()  # number of channels
    f.readline()  # number of dye seqs
    weightmap = {}
    while (True):
        line = f.readline()
        if (line == ""):
            break
        ltabs = line.split("\t")
        weightmap[int(ltabs[2])] = int(ltabs[1])
    f.close()
    f = open(directory + predictions_file, "r")
    csv = f.readlines()
    csv = csv[1:]
    predictions = [0] * len(csv)
    for i in range(0, len(csv)):
        cells = csv[i].split(",")
        predictions[i] = (int(cells[1]), float(cells[2]))
    f.close()
    precision, recall = compute_pr_curve(predictions, true_y, weightmap)
    ax.plot(recall, precision, '-', label = label, linewidth = 2)

def compute_pr_curve(predictions, ground_truth, weightmap):
    amt_correct_and_score = [(0.0, 0.0)] * len(predictions)

    correct = 0.0
    incorrect = 0.0
    for i in range(len(ground_truth)):
        if predictions[i][0] == ground_truth[i]:
            amt_correct = 1.0 / weightmap[predictions[i][0]]
            amt_correct_and_score[i] = (amt_correct, predictions[i][1])
            correct += amt_correct
            incorrect += 1.0 - amt_correct
        else:
            amt_correct_and_score[i] = (0.0, predictions[i][1])
            incorrect += 1.0
#    print('precision at highest recall: ' + str(correct) + '/'
#          + str(correct + incorrect) + ', or '
#          + str(correct / (correct + incorrect)))

    correct = 0
    incorrect = 0
    total = len(predictions)
    precision = []
    recall = []
    amt_correct_and_score.sort(key = lambda entry: entry[1], reverse = True)
    for i in range(len(amt_correct_and_score)):
        entry = amt_correct_and_score[i]
        next_entry = 0
        if i + 1 < len(amt_correct_and_score):
            next_entry = amt_correct_and_score[i + 1]
        else:
            next_entry = (0.0, -1)
        correct += entry[0]
        incorrect += 1.0 - entry[0]
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
