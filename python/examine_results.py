# -*- coding: utf-8 -*-
"""
@author: Matthew Beauregard Smith (UT Austin)
"""

def examine_results(predictions_file, radiometries_file, directory):
    f = open(directory + predictions_file, "r")
    preds_csv = f.readlines()
    f.close()
    preds_csv = preds_csv[1:]
    preds = [0] * len(preds_csv)
    for i in range(0, len(preds_csv)):
        cells = preds_csv[i].split(",")
        preds[i] = (int(cells[1]), float(cells[2]))
    f = open(directory + radiometries_file, "r")
    rads_tsv = f.readlines()
    f.close()
    rads_tsv = rads_tsv[3:]
    for i in range(0, len(preds)):
        cells = rads_tsv[i].split("\t")
        preds[i] = (preds[i][0], preds[i][1], cells[0:18:3])
    preds.sort(key = lambda entry: entry[1], reverse = True)
    for i in range(100):
        print(preds[i])
