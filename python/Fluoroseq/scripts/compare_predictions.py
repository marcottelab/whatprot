# -*- coding: utf-8 -*-
"""
@author: Matthew Beauregard Smith (UT Austin)
"""

# This file is meant to be used with a MODIFIED version of the nearest
# neighbors code. It does not work on any .csv files.

from statistics import mean

HMM_FILE = "C:/Users/Matthew/ICES/MarcotteLab/data/classification/control_15_proteins/best_hmm_results.csv"
NN_FILE_1ST = "C:/Users/Matthew/ICES/MarcotteLab/data/classification/control_15_proteins/ann_1st.csv"
NN_FILE_2ND = "C:/Users/Matthew/ICES/MarcotteLab/data/classification/control_15_proteins/ann_2nd.csv"
NN_FILE_3RD = "C:/Users/Matthew/ICES/MarcotteLab/data/classification/control_15_proteins/ann_3rd.csv"
NN_FILE_4TH = "C:/Users/Matthew/ICES/MarcotteLab/data/classification/control_15_proteins/ann_4th.csv"
NN_FILE_5TH = "C:/Users/Matthew/ICES/MarcotteLab/data/classification/control_15_proteins/ann_5th.csv"

f = open(HMM_FILE, 'r')
hmm_csv = f.readlines()
hmm_csv = hmm_csv[1:]
hmm_ids = [0] * len(hmm_csv)
for i in range(len(hmm_csv)):
    hmm_ids[i] = int(hmm_csv[i].split(",")[1])
f.close()

f = open(NN_FILE_1ST, 'r')
nn_csv_1st = f.readlines()
nn_csv_1st = nn_csv_1st[1:]
nn_ids_1st = [0] * len(nn_csv_1st)
for i in range(len(nn_csv_1st)):
    nn_ids_1st[i] = int(nn_csv_1st[i].split(",")[1])
f.close()

f = open(NN_FILE_2ND, 'r')
nn_csv_2nd = f.readlines()
nn_csv_2nd = nn_csv_2nd[1:]
nn_ids_2nd = [0] * len(nn_csv_2nd)
for i in range(len(nn_csv_2nd)):
    nn_ids_2nd[i] = int(nn_csv_2nd[i].split(",")[1])
f.close()

f = open(NN_FILE_3RD, 'r')
nn_csv_3rd = f.readlines()
nn_csv_3rd = nn_csv_3rd[1:]
nn_ids_3rd = [0] * len(nn_csv_3rd)
for i in range(len(nn_csv_3rd)):
    nn_ids_3rd[i] = int(nn_csv_3rd[i].split(",")[1])
f.close()

f = open(NN_FILE_4TH, 'r')
nn_csv_4th = f.readlines()
nn_csv_4th = nn_csv_4th[1:]
nn_ids_4th = [0] * len(nn_csv_4th)
for i in range(len(nn_csv_4th)):
    nn_ids_4th[i] = int(nn_csv_4th[i].split(",")[1])
f.close()

f = open(NN_FILE_5TH, 'r')
nn_csv_5th = f.readlines()
nn_csv_5th = nn_csv_5th[1:]
nn_ids_5th = [0] * len(nn_csv_5th)
for i in range(len(nn_csv_5th)):
    nn_ids_5th[i] = int(nn_csv_5th[i].split(",")[1])
f.close()

num_matches_1st = 0
num_matches_2nd = 0
num_matches_3rd = 0
num_matches_4th = 0
num_matches_5th = 0
for i in range(len(hmm_ids)):
    num_matches_1st += int(hmm_ids[i] == nn_ids_1st[i])
    num_matches_2nd += int(hmm_ids[i] == nn_ids_2nd[i])
    num_matches_3rd += int(hmm_ids[i] == nn_ids_3rd[i])
    num_matches_4th += int(hmm_ids[i] == nn_ids_4th[i])
    num_matches_5th += int(hmm_ids[i] == nn_ids_5th[i])
num_matches_all = num_matches_1st + num_matches_2nd + num_matches_3rd + num_matches_4th + num_matches_5th

print("total count: " + str(len(hmm_ids)))
print("")
print("1st, matches: " + str(num_matches_1st) + ",    rate: " + str(num_matches_1st / len(hmm_ids)))
print("2nd, matches: " + str(num_matches_2nd) + ",    rate: " + str(num_matches_2nd / len(hmm_ids)))
print("3rd, matches: " + str(num_matches_3rd) + ",    rate: " + str(num_matches_3rd / len(hmm_ids)))
print("4th, matches: " + str(num_matches_4th) + ",    rate: " + str(num_matches_4th / len(hmm_ids)))
print("5th, matches: " + str(num_matches_5th) + ",    rate: " + str(num_matches_5th / len(hmm_ids)))
print("")
print("all, matches: " + str(num_matches_all) + ",    rate: " + str(num_matches_all / len(hmm_ids)))
