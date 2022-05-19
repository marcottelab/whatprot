# -*- coding: utf-8 -*-
"""
@author: Matthew B. Smith (UT Austin)
"""

def fix_ids(true_id_file, full_peps_file, lim_peps_file, fixed_true_id_file):
    f = open(full_peps_file, 'r')
    prot2pepid = {}
    curprot = -1
    pepid = 0
    while True:
        line = f.readline()
        if len(line) <= 0:
            break
        protid = int(line.split(',')[1])
        if protid != curprot:
            prot2pepid[protid] = pepid
            curprot = protid
        pepid += 1
    f.close()
    f = open(lim_peps_file, 'r')
    badid2fixedid = {}
    curprot = -1
    protoffset = 0
    badid = 0
    while True:
        line = f.readline()
        if len(line) <= 0:
            break
        protid = int(line.split(',')[1])
        firstpep = prot2pepid[protid]
        if protid != curprot:
            protoffset = 0
            curprot = protid
        else:
            protoffset += 1
        badid2fixedid[badid] = firstpep + protoffset
        badid += 1
    f.close()
    f = open(true_id_file, 'r')
    g = open(fixed_true_id_file, 'w')
    g.write(f.readline())  # Transfer number of rows to new file
    while True:
        line = f.readline()
        if len(line) <= 0:
            break
        badid = int(line)
        fixedid = badid2fixedid[badid]
        g.write(str(fixedid) + '\n')
    f.close()
    g.close()
