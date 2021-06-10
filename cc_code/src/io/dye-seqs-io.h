/******************************************************************************\
* Author: Matthew Beauregard Smith                                             *
* Affiliation: The University of Texas at Austin                               *
* Department: Oden Institute and Institute for Cellular and Molecular Biology  *
* PI: Edward Marcotte                                                          *
* Project: Protein Fluorosequencing                                            *
\******************************************************************************/

#ifndef WHATPROT_IO_DYE_SEQS_IO_H
#define WHATPROT_IO_DYE_SEQS_IO_H

// Standard C++ library headers:
#include <string>
#include <vector>

// Local project headers:
#include "common/dye-seq.h"
#include "common/sourced-data.h"

namespace whatprot {

void read_dye_seqs(
        const std::string& filename,
        unsigned int* num_channels,
        unsigned int* total_num_dye_seqs,
        std::vector<SourcedData<DyeSeq, SourceCount<int>>>* dye_seqs);

void read_dye_seqs_raw(const std::string& filename,
                       unsigned int* num_channels,
                       unsigned int* num_dye_seqs,
                       unsigned int** dye_string_lengths,
                       char*** dye_strings,
                       unsigned int** dye_seqs_num_peptides,
                       int** dye_seqs_ids);

void convert_dye_seqs_from_raw(
        unsigned int num_channels,
        unsigned int num_dye_seqs,
        unsigned int* dye_string_lengths,
        char** dye_strings,
        unsigned int* dye_seqs_num_peptides,
        int* dye_seqs_ids,
        std::vector<SourcedData<DyeSeq, SourceCount<int>>>* dye_seqs);

}  // namespace whatprot

#endif  // WHATPROT_IO_DYE_SEQS_IO_H
