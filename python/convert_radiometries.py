# -*- coding: utf-8 -*-
"""
@author: Matthew Beauregard Smith (UT Austin)
"""

import numpy as np

# num_channels: number of channels (or colors) of fluorophore in use.
# num_mocks: number of Edman mocks before sequencing truly begins. These cycles
#            are discarded.
# num_cycles: number of images taken during sequencing. This is one higher than
#             the comparable value from Erisyon, because they call the image
#             taken before the first Edman a 'pre,' and every run will have
#             1 pre and x Edmans. So set num_cycles to x+1.
# radmat_file: input file from Erisyon, in either .tsv or .npy format.
# output_file: where you want your data to go. The filename should end in .tsv,
#              since it is produced in .tsv (tab separated values) format.
def convert_radiometries(num_channels, num_mocks, num_cycles, radmat_file, output_file, ranges_by_channel=None):
    # radmats may come in two different formats, a .npy (numpy) file or a .tsv
    # (tab separated values) file. We can handle either but we need to know the
    # extension.
    radmat = None
    radmat_filetype = radmat_file.split('.')[-1]
    if radmat_filetype == 'npy':
        radmat = np.load(radmat_file)
    elif radmat_filetype == 'tsv':
        radmat = np.genfromtxt(radmat_file, delimiter='\t', dtype=float)[:,1:]
        radmat = np.reshape(radmat, (radmat.shape[0], num_channels, num_mocks + num_cycles))
    else:
        # This is a problem.
        return
    
    print('Read ' + str(radmat.shape[0]) + ' radiometries.\n')

    # Filter out user defined ranges
    if ranges_by_channel != None:
        row_list = []
        for i in range(radmat.shape[0]):
            row_ok = True
            for c in range(num_channels):
                channel_ok = True
                for j in range(num_mocks + num_cycles):
                    x = radmat[i, c, j]
                    value_ok = False
                    for k in range(len(ranges_by_channel[c])):
                        if ranges_by_channel[c][k][0] < x and x < ranges_by_channel[c][k][1]:
                            value_ok = True
                            break
                    if not value_ok:
                        channel_ok = False
                        break
                if not channel_ok:
                    row_ok = False
                    break
            if row_ok:
                row_list += [i]
        radmat = radmat[row_list,:,:]
        print('After filtering, ' + str(radmat.shape[0]) + ' radiometries remain.\n')

    radmat = np.transpose(radmat, (0, 2, 1))

    # Fix intensities. If the intensities are systematically too large, the
    # probabilities will be too low and we will get underflow when the values
    # are fed into whatprot. Switching into log-probability format is not a
    # very viable option, because a significant amount of addition is needed,
    # and addition is very slow in log-probability format. However the intensity
    # values are of arbitrary units anyways. We have found that when working
    # with Erisyon data, dividing by 15000 is an appropriate adjustment, but
    # this number was rather arbitrarily chosen.
    radmat = radmat / 15000

    f = open(output_file, 'w')
    # The output file starts with three lines of metadata useful for
    # understanding the specifics of this file. These are in a standardized
    # ordering.
    f.write(str(num_cycles) + "\n")
    f.write(str(num_channels) + "\n")
    f.write(str(radmat.shape[0]) + "\n")
    for rad in radmat:
        # We omit all mock cycles. Notice also that one row of output includes
        # both the cycle dimension AND the channels dimension.
        for i in range(num_mocks, num_mocks + num_cycles):
            for j in range(num_channels):
                # tab between all entries but not before first or after last in
                # row.
                if i > 0 or j > 0:
                    f.write("\t")
                f.write(str(rad[i, j]))
        f.write("\n")
    f.close()
