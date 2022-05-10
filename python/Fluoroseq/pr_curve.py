# -*- coding: utf-8 -*-
"""
@author: Matthew Beauregard Smith (UT Austin)
"""

import matplotlib.pyplot as plt
from collections import defaultdict
from numpy import load

plt.rcParams["font.family"] = "Calibri"

def plot_pr_curve(predictions_file,
                  true_y_file,
                  dye_seqs_file,
                  full_peps_file = None,
                  lim_peps_file = None,
                  directory = ""):
    plot_pr_curves([predictions_file],
                   true_y_file,
                   dye_seqs_file,
                   full_peps_file,
                   lim_peps_file,
                   directory)

def plot_pr_curves(predictions_files,
                   true_y_file,
                   dye_seqs_file,
                   full_peps_file = None,
                   lim_peps_file = None,
                   directory = "",
                   labels = []):
    fig = 0
    axs = 0
    if full_peps_file == None:
        fig, ax = plt.subplots(1)
        axs = [ax]
    else:
        fig, axs = plt.subplots(3)
    cmloc = [x / len(predictions_files) for x in range(len(predictions_files))]
    for ax in axs:
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
        cur_full_peps_file = ''
        if full_peps_file != None:
            if isinstance(full_peps_file, str):
                cur_full_peps_file = full_peps_file
            else:
                cur_full_peps_file = full_peps_file[i]
        cur_lim_peps_file = ''
        if lim_peps_file != None:
            if isinstance(lim_peps_file, str):
                cur_lim_peps_file = lim_peps_file
            else:
                cur_lim_peps_file = lim_peps_file[i]
        plot_pr_curve_noshow(axs[0],
                             predictions_file,
                             cur_true_y_file,
                             cur_dye_seqs_file,
                             directory = directory,
                             label = label)
        if full_peps_file != None:
            plot_aggregate_pr_curves_noshow(axs[1],
                                            axs[2],
                                            predictions_file,
                                            cur_dye_seqs_file,
                                            cur_full_peps_file,
                                            cur_lim_peps_file,
                                            directory = directory,
                                            label = label)
    for ax in axs:
        ax.set_xlim(left = 0)
        ax.set_ylim(bottom = 0)
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        ax.spines["bottom"].set_linewidth(1.5)
        ax.spines["left"].set_linewidth(1.5)
        ax.xaxis.set_tick_params(width=1.5)
        ax.yaxis.set_tick_params(width=1.5)
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
        id_or_ids = eval(ltabs[2])
        dseq_id = -1
        if type(id_or_ids) is int:
            dseq_id = id_or_ids
        else:
            dseq_id = id_or_ids[0]
        weightmap[dseq_id] = int(ltabs[1])
    f.close()
    f = open(directory + predictions_file, "r")
    csv = f.readlines()
    csv = csv[1:]
    predictions = [0] * len(csv)
    for i in range(0, len(csv)):
        cells = csv[i].split(",")
        predictions[i] = (int(cells[1]), float(cells[2]))
    f.close()
    precision, recall = compute_read_pr_curve(predictions, true_y, weightmap)
    ax.plot(recall, precision, '-', label = label, linewidth = 2)

def plot_aggregate_pr_curves_noshow(axpep,
                                    axprot,
                                    predictions_file,
                                    dye_seqs_file,
                                    full_peps_file,
                                    lim_peps_file,
                                    directory = "",
                                    title = "Precision / Recall",
                                    label = None):
    f = open(directory + full_peps_file, 'r')
    csv = f.readlines()
    pep2prot = {}
    prot2peps = defaultdict(list)
    for i in range(len(csv)):
        cells = csv[i].split(',')
        protid = int(cells[1])
        pep2prot[i] = protid
        prot2peps[protid] += [i]
    f.close()
    f = open(directory + lim_peps_file, 'r')
    csv = f.readlines()
    prots = defaultdict(bool)  # defaults to False
    for i in range(len(csv)):
        cells = csv[i].split(',')
        prots[int(cells[1])] = True
    f.close()
    peps = defaultdict(bool)  # defaults to False
    for protid, _ in prots.items():
        for pepid in prot2peps[protid]:
            peps[pepid] = True
    f = open(directory + dye_seqs_file, 'r')
    tsv = f.readlines()
    tsv = tsv[2:]
    dseqs2peps = defaultdict(list)  # defaults to []
    for i in range(len(tsv)):
        cells = tsv[i].split('\t')
        # cells[0] is the actual sequence; we don't care.
        # dseq id is i
        # cells[1] is the number of peps; it's redundant since cells[2] is list
        dseq_peps = eval(cells[2])
        # 0th pep is used as dseq id as shortcut for other purposes.... ew.
        dseqs2peps[dseq_peps[0]] = dseq_peps
    f.close()
    f = open(directory + predictions_file, 'r')
    csv = f.readlines()
    csv = csv[1:]
    pep2score = defaultdict(float)  # defaults to 0.0
    for i in range(len(csv)):
        cells = csv[i].split(',')
        # cells[0] is radiometry id, same as i
        dseq_id = int(cells[1])
        pep_score = float(cells[2])  # Note: score is already corrected for hits count
        for pep_id in dseqs2peps[dseq_id]:
            pep2score[pep_id] = max(pep2score[pep_id], pep_score)
    f.close()
    pep_precision, pep_recall = compute_aggregate_pr_curve(pep2score, peps)
    axpep.plot(pep_recall, pep_precision, '-', label = label, linewidth = 2)
    prot2score = defaultdict(float)  # defaults to 0.0
    for pep, pepscore in pep2score.items():
        prot = pep2prot[pep]
        prot2score[prot] = max(prot2score[prot], pepscore)
    prot_precision, prot_recall = compute_aggregate_pr_curve(prot2score, prots)
    axprot.plot(prot_recall, prot_precision, '-', label = label, linewidth = 2)

def compute_aggregate_pr_curve(predictions, ground_truth):
    len_ground_truth = len(ground_truth)  # must do this before access of elems with False messes up result.
    amt_correct_and_score = [(0.0, 0.0)] * len(predictions)
    i = 0
    for pred_id, pred_score in predictions.items():
        if ground_truth[pred_id]:
            amt_correct_and_score[i] = (1.0, pred_score)
        else:
            amt_correct_and_score[i] = (0.0, pred_score)
        i += 1
    return compute_pr_curve_helper(amt_correct_and_score, len_ground_truth)

def compute_read_pr_curve(predictions, ground_truth, weightmap):
    amt_correct_and_score = [(0.0, 0.0)] * len(predictions)
    for i in range(len(ground_truth)):
        if predictions[i][0] == ground_truth[i]:
            amt_correct = 1.0 / weightmap[predictions[i][0]]
            amt_correct_and_score[i] = (amt_correct, predictions[i][1])
        else:
            amt_correct_and_score[i] = (0.0, predictions[i][1])
    return compute_pr_curve_helper(amt_correct_and_score, len(amt_correct_and_score))

def compute_pr_curve_helper(amt_correct_and_score, total):
    correct = 0.0
    incorrect = 0.0
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
            precision.append([correct / (correct + incorrect)])
            recall.append(correct / total)
    for i in range(len(recall)):
        precision[i] = sum(precision[i]) / len(precision[i])
    return (precision, recall)
