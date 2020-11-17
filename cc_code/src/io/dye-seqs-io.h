/******************************************************************************\
* Author: Matthew Beauregard Smith                                             *
* Affiliation: The University of Texas at Austin                               *
* Department: Oden Institute and Institute for Cellular and Molecular Biology  *
* PI: Edward Marcotte                                                          *
* Project: Protein Fluorosequencing                                            *
\******************************************************************************/

// For MPI version, define compiler macro USE_MPI when building.

#ifndef FLUOROSEQ_IO_DYE_SEQS_IO_H
#define FLUOROSEQ_IO_DYE_SEQS_IO_H

// Standard C++ library headers:
#include <string>
#include <vector>

// MPI header:
#ifdef USE_MPI
#include <mpi.h>
#endif  // USE_MPI

// Local project headers:
#include "common/dye-seq.h"
#include "common/sourced-data.h"

namespace fluoroseq {

enum MpiReadMode { BROADCAST, SCATTER };

void read_dye_seqs(const std::string& filename,
                   int* num_channels,
                   int* total_num_dye_seqs,
                   std::vector<SourcedData<DyeSeq, SourceCount<int>>>* dye_seqs,
                   MpiReadMode mode);

void read_dye_seqs_raw(const std::string& filename,
                       int* num_channels,
                       int* num_dye_seqs,
                       int** dye_string_lengths,
                       char*** dye_strings,
                       int** dye_seqs_num_peptides,
                       int** dye_seqs_ids);

#ifdef USE_MPI

void broadcast_dye_seqs(int* num_channels,
                        int* num_dye_seqs,
                        int** dye_string_lengths,
                        char*** dye_strings,
                        int** dye_seqs_num_peptides,
                        int** dye_seqs_ids);

void scatter_dye_seqs(int* num_channels,
                      int* total_num_dye_seqs,
                      int* num_dye_seqs,
                      int** dye_string_lengths,
                      char*** dye_strings,
                      int** dye_seqs_num_peptides,
                      int** dye_seqs_ids);

#endif  // USE_MPI

void convert_dye_seqs_from_raw(
        int num_channels,
        int num_dye_seqs,
        int* dye_string_lengths,
        char** dye_strings,
        int* dye_seqs_num_peptides,
        int* dye_seqs_ids,
        std::vector<SourcedData<DyeSeq, SourceCount<int>>>* dye_seqs);

}  // namespace fluoroseq

#endif  // FLUOROSEQ_IO_DYE_SEQS_IO_H
