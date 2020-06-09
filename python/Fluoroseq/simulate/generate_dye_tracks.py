# -*- coding: utf-8 -*-
"""
@author: Matthew B. Smith (UT Austin)
"""

from collections import defaultdict
from common.dye_track import DyeTrack
from numpy import array
from simulate.simulate_sequencing import simulate_sequencing

def generate_dye_tracks(error_model,
                        n_channels,
                        n_timesteps,
                        dye_seqs,
                        n_samples_per_dye_seq):
    em = error_model.copy()
    em.sigma = 0
    em.bg_lambda = 0
    dtdict = defaultdict(list)
    for i in range(len(dye_seqs)):
        for j in range(n_samples_per_dye_seq):
            radiometry = simulate_sequencing(em,
                                             n_channels,
                                             n_timesteps,
                                             dye_seqs[i])
            dye_track = DyeTrack(array(radiometry.intensities, dtype=int),
                                 radiometry.src_dye_seq)
            dtdict[str(dye_track)] += [dye_track]
    dye_tracks = [0] * len(dtdict)
    index = 0
    for (_, value) in dtdict.items():
        dye_track = DyeTrack(value[0].intensities, defaultdict(int))
        for dt in value:
            dye_track.src_dye_seq[dt.src_dye_seq] += 1
        dye_tracks[index] = dye_track
        index += 1
    return dye_tracks