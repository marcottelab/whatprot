# -*- coding: utf-8 -*-
"""
@author: Matthew Beauregard Smith (UT Austin)
"""

from common.error_model import ErrorModel
from common.peptide import Peptide
from simulate.label_peptides import label_peptides
from simulate.generate_dye_tracks import generate_dye_tracks

NUM_PEPTIDES = 705
NUM_CHANNELS = 3
NUM_CYCLES = 16
NUM_SAMPLES_PER_PEPTIDE = 1000
LABEL_SET = ['DE','Y','C']
PEPTIDE_FILE = 'C:/Users/Matthew/ICES/MarcotteLab/data/classification/control_15_proteins/peps.csv'
OUTPUT_FILE = 'C:/Users/Matthew/ICES/MarcotteLab/data/classification/control_15_proteins/dye_tracks_small.tsv'

f = open(PEPTIDE_FILE, 'r')
f.readline()  # header
f.readline()  # Zack's null line
line = f.readline()
peptides = [0] * NUM_PEPTIDES
i = 0
while line != '\n' and line != '':
    items = line.split(",")
    pep_id = items[0]
    pep_str = items[-1]
    peptides[i] = Peptide(pep_str, pep_id=pep_id)
    line = f.readline()
    i += 1
f.close()
dye_seqs = label_peptides(peptides, LABEL_SET)

em = ErrorModel(.94, .05, .07, 1, .16, 0, .05, True)
dye_tracks = generate_dye_tracks(em,
                                 NUM_CHANNELS,
                                 NUM_CYCLES,
                                 dye_seqs,
                                 NUM_SAMPLES_PER_PEPTIDE)

f = open(OUTPUT_FILE, 'w')
f.write(str(NUM_CYCLES) + "\n")
f.write(str(NUM_CHANNELS) + "\n")
f.write(str(len(dye_tracks)) + "\n")
for dye_track in dye_tracks:
    f.write(str(dye_track.feature_vector()[0]))
    for i in range(1, NUM_CHANNELS * NUM_CYCLES):
        f.write("\t" + str(dye_track.feature_vector()[i]))
    f.write("\t" + str(len(dye_track.src_dye_seq)))
    for (key, value) in dye_track.src_dye_seq.items():
        f.write("\t" + str(key.src_peptides[0].pep_id))
        f.write("\t" + str(len(key.src_peptides)))
        f.write("\t" + str(value))
    f.write("\n")
f.close()


