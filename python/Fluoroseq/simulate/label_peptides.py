# -*- coding: utf-8 -*-
"""
@author: Matthew B. Smith (UT Austin)
"""

from collections import defaultdict
from common.dye_seq import DyeSeq

def label_peptides(peptides, labels_by_channel):
    str_dye_seqs = defaultdict(lambda: [])
    for peptide in peptides:
        list_dye_seq = ['.'] * len(peptide.seq)
        for i in range(len(peptide.seq)):
            for ch in range(len(labels_by_channel)):
                for aa in labels_by_channel[ch]:
                    if aa == peptide.seq[i]:
                        list_dye_seq[i] = str(ch)
        while len(list_dye_seq) > 0 and list_dye_seq[-1] == '.':
            list_dye_seq.pop()
        str_dye_seq = ''.join(list_dye_seq)
        if str_dye_seq != '':
            str_dye_seqs[str_dye_seq] += [peptide]
    dye_seqs = [0] * len(str_dye_seqs)
    index = 0
    for str_dye_seq in str_dye_seqs:
        dye_seqs[index] = DyeSeq(len(labels_by_channel), str_dye_seq, index,
                str_dye_seqs[str_dye_seq])
        index += 1
    return dye_seqs
#
#def label_peptide(peptide, labels_by_channel, pep_id):
#    list_dye_seq = ['.'] * len(peptide.seq)
#    for i in range(len(peptide.seq)):
#        for ch in range(len(labels_by_channel)):
#            for aa in labels_by_channel[ch]:
#                if aa == peptide.seq[i]:
#                    list_dye_seq[i] = str(ch)
#    while len(list_dye_seq) > 0 and list_dye_seq[-1] == '.':
#        list_dye_seq.pop()
#    str_dye_seq = ''.join(list_dye_seq)
#    dye_seq = DyeSeq(len(labels_by_channel), str_dye_seq, pep_id, [peptide])
#    return dye_seq
