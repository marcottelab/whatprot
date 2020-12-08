# -*- coding: utf-8 -*-
"""
@author: Matthew Beauregard Smith (UT Austin)
"""

def trim_lines(filein, fileout, directory, mod, offset = 0):
    fin = open(directory + filein, 'r')
    fout = open(directory + fileout, 'w')
    i = 1
    while (i <= offset):
        fin.readline()
        i += 1
    i = 1
    while (True):
        line = fin.readline()
        if (line == ""):
            break
        if (i % mod == 1):
            fout.write(line)
        i += 1
    fin.close()
    fout.close()
