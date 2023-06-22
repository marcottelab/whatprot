# -*- coding: utf-8 -*-
"""
@author: Matthew Beauregard Smith (UT Austin)
"""

import matplotlib.pyplot as plt
import matplotlib.ticker as mtick
import numpy as np
import math
# (bw_bootstrap_file, bw_best, powell_bootstrap_file, powell_best)
def bootstrap_fit_boxplot(files, path, bests, doubled=True, onefluor=False, target_values=None):
    csvs = [0] * len(files)
    if onefluor:
        for i in range(len(files)):
            if i % 2 == 0 or not doubled:
                csvs[i] = np.loadtxt(path + files[i], skiprows=1, usecols=(0,2,4), delimiter=',', dtype=float)
            else:
                csvs[i] = np.loadtxt(path + files[i], skiprows=1, usecols=(4,3,5), delimiter=',', dtype=float)
    else:
        for i in range(len(files)):
            if i % 2 == 0 or not doubled:
                csvs[i] = np.loadtxt(path + files[i], skiprows=1, usecols=(0,1,2,3,4,5), delimiter=',', dtype=float)
            else:
                csvs[i] = np.loadtxt(path + files[i], skiprows=1, usecols=(6,5,3,4,7,8), delimiter=',', dtype=float)
    csvs = np.stack(csvs)
    numcols = 6
    if onefluor == True:
        numcols = 3
    for i in range(numcols):
        plt.subplot(1,numcols,i+1)
        bp = plt.boxplot(csvs[:,:,i].transpose(), patch_artist=True, widths=0.75, medianprops=dict(linewidth=2), flierprops=dict(marker='.', markersize=2, zorder=11))
        num_exp = len(files)
        colors = []
        lcolors = []
        wcolors = []
        mcolors = []
        linestyles = []
        hatches = []
        if doubled:
            num_exp //= 2
            colors = ['dodgerblue','lightcoral'] * num_exp
            lcolors = ['darkblue', 'darkred'] * num_exp
            wcolors = ['darkblue', 'darkblue', 'darkred', 'darkred'] * num_exp
            linestyles = ['-','-','-','-','--','--','--','--',':',':',':',':','-.','-.','-.','-.']
            hatches = ['', '', '////', '////', '....', '....', 'OO', 'OO']
        else:
            colors = ['dodgerblue'] * num_exp
            lcolors = ['darkblue'] * num_exp
            wcolors = ['darkblue', 'darkblue'] * num_exp
            linestyles = ['-','-','--','--',':',':','-.','-.']
            hatches = ['', '////', '....', 'OO']
        # Must set color before setting facecolor (yes, really)
        for patch, color in zip(bp['boxes'], lcolors):
            patch.set_color(color) # box border + hatch color (must be the same)
        for patch, color in zip(bp['boxes'], colors):
            patch.set_facecolor(color) # box fill color
        for patch, color in zip(bp['fliers'], lcolors):
            patch.set_markeredgecolor(color) # flier color (no fill)
        for patch, color, linestyle in zip(bp['whiskers'], wcolors, linestyles):
            patch.set_color(color)
            patch.set_linestyle(linestyle)
        for patch, color in zip(bp['caps'], wcolors):
            patch.set_color(color)
        for patch, hatch in zip(bp['boxes'], hatches):
            patch.set_hatch(hatch)
        for patch, color in zip(bp['medians'], lcolors):
            patch.set_color(color)
        for j in range(len(bests)):
            color = 'darkblue'
            if doubled and j % 2 == 1:
                color = 'darkred'
            plt.plot(j + 0.75, bests[j][i], color='white', markeredgecolor=color, markeredgewidth='1', marker='>', zorder=10)
        if target_values:
            plt.plot((0,1+len(files)),(target_values[i],target_values[i]), linestyle='--', color='black')
        plt.gca().yaxis.set_major_formatter(mtick.PercentFormatter(1))
        plt.ylim((math.floor(min(csvs[:,:,i].flatten())*100)/100, math.ceil(max(csvs[:,:,i].flatten())*100)/100))
        plt.tick_params(axis='x',which='both',bottom=False,top=False,labelbottom=False)
    plt.tight_layout(pad=0)
