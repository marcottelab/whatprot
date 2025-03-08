# -*- coding: utf-8 -*-
"""
@author: Matthew Beauregard Smith (UT Austin)
"""

from common.peptide import Peptide
from simulate.label_peptides import label_peptides

def dye_seqs_from_peptides(peptide_file, label_set, dye_seqs_file, peps_prots_file, amts = None, mode='first'):
    f = open(peptide_file, 'r')
    # f.readline()  # header
    # f.readline()  # Zack's null line
    line = "placeholder"
    peptides = []
    proteins = set()
    pep_id = 0
    pep2prot = {}
    while line != '\n' and line != '':
        line = f.readline()[0 : -1]
        pattern = line.split(",")[0]
        peptides += [Peptide(pattern, pep_id=pep_id)]
        prot_id = int(line.split(",")[1])
        proteins.add(prot_id)
        pep2prot[pep_id] = prot_id
        pep_id += 1
    f.close()
    dye_seqs = label_peptides(peptides, label_set)
    
    f = open(dye_seqs_file, 'w')
    f.write(str(len(label_set)) + "\n")  # num channels
    f.write(str(len(dye_seqs)) + "\n")
    for dye_seq in dye_seqs:
        dye_seq.dye_seq.reverse()
        f.write("".join(dye_seq.dye_seq) + "\t")
        f.write(str(len(dye_seq.src_peptides)) + "\t")
        if mode == 'first':
            f.write(str(dye_seq.src_peptides[0].pep_id) + "\n")
        elif mode == 'all':
            f.write(str([x.pep_id for x in dye_seq.src_peptides]) + "\n")
        else:
            print("bad-input, invalid mode: " + mode)
    f.close()

    f = open(peps_prots_file, 'w')
    proteins = list(proteins).sort()
    f.write('protein,' + ','.join(proteins) + '\n')
    f.write('amt,' + ','.join(amts) + '\n')
    f.write('peptide\n')
    for dye_seq in dye_seqs:
        f.write()

