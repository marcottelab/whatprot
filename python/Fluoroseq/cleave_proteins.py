# -*- coding: utf-8 -*-
"""
@author: Matthew B. Smith (UT Austin)
"""

from random import sample

def cleave_proteins(fasta, fpeptides, protease, n = -1):
    fpro = open(fasta, "r")
    fpro.readline()  # skip first '>' line for convenience.
    proteins = []
    while True:
        protein = ""
        line = ""
        while True:
            line = fpro.readline()[0 : -1]
            if (not line):
                break
            if (line[0] == '>'):
                break
            protein += line
        if (not line):
            break
        proteins += [protein]
    fpro.close()
    npros = 0
    if (n == -1):
        npros = proteins
    else:
        npros = sample(proteins, n)
    print(len(npros))
    fpep = open(fpeptides, "w")
    # Conventionally fasta files have proteins from N-terminus to C-terminus
    if (protease == "trypsin"):
        # Now we trypsinize, cutting the C-terminal side of lysine (K) and
        # arginine (R) residues, but only if they are not followed by Proline
        # (P).
        for protein in npros:
            lastcut = 0
            for i in range(len(protein) - 1):
                if ((protein[i] == 'K' or protein[i] == 'R') and (protein[i + 1] != 'P')):
                    peptide = protein[lastcut : i + 1]
                    lastcut = i + 1
                    fpep.write(peptide + "\n")
            peptide = protein[lastcut:]
            fpep.write(peptide + "\n")
    elif (protease == "cyanogen bromide"):
        # Now we use cyanogen bromide, cutting the C-terminal side of methionine
        # (M).
        for protein in npros:
            lastcut = 0
            for i in range(len(protein) - 1):
                if ((protein[i] == 'M')):
                    peptide = protein[lastcut : i + 1]
                    lastcut = i + 1
                    fpep.write(peptide + "\n")
            peptide = protein[lastcut:]
            fpep.write(peptide + "\n")
    fpep.close()
