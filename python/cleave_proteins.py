# -*- coding: utf-8 -*-
"""
@author: Matthew B. Smith (UT Austin)
"""

from random import sample

def cleave_proteins(fasta, fpeptides, protease, n = -1):
    fpro = open(fasta, "r")
    fpro.readline()  # skip first '>' line for convenience.
    proteins = []
    pid = 0
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
        proteins += [(pid, protein)]
        if (not line):
            break
        pid += 1
    fpro.close()
    npros = 0
    if (n == -1):
        npros = proteins
    else:
        npros = sample(proteins, n)
    print(len(npros))
    fpep = open(fpeptides, "w")
    # Conventionally fasta files have proteins from N-terminus to C-terminus
    for (pid, protein) in npros:
        lastcut = 0
        for i in range(len(protein) - 1):
            cut = False
            if (protease == "trypsin"):
                # Here we trypsinize, cutting the C-terminal side of lysine (K)
                # and arginine (R) residues, but only if they are not followed
                # by Proline (P).
                if ((protein[i] == 'K' or protein[i] == 'R') and (protein[i + 1] != 'P')):
                    cut = True
            elif (protease == "cyanogen bromide"):
                # Now we use cyanogen bromide, cutting the C-terminal side of
                # methionine (M).
                if ((protein[i] == 'M')):
                    cut = True
            elif (protease == "EndoPRO"):
                # Here we use EndoPRO, cutting the C-terminal side of proline
                # (P) and alanine (A).
                if ((protein[i] == 'P' or protein[i] == 'A')):
                    cut = True
            else:
                print('error, invalid protease: ' + protease)
            if (cut == True):
                peptide = protein[lastcut : i + 1]
                lastcut = i + 1
                fpep.write(peptide + "," + str(pid) + "\n")
        peptide = protein[lastcut:]
        fpep.write(peptide + "," + str(pid) + "\n")
    fpep.close()
