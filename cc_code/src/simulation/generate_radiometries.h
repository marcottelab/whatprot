/******************************************************************************\
* Author: Matthew Beauregard Smith                                             *
* Affiliation: The University of Texas at Austin                               *
* Department: Oden Institute and Institute for Cellular and Molecular Biology  *
* PI: Edward Marcotte                                                          *
* Project: Protein Fluorosequencing                                            *
\******************************************************************************/

#ifndef FLUOROSEQ_SIMULATION_GENERATE_RADIOMETRIES_H
#define FLUOROSEQ_SIMULATION_GENERATE_RADIOMETRIES_H

// Standard C++ library headers:
#include <random>
#include <vector>

// Local project headers:
#include "common/dye_seq.h"
#include "common/error_model.h"
#include "common/radiometry.h"
#include "common/sourced_data.h"

namespace fluoroseq {

void generate_radiometries(
        const ErrorModel& error_model,
        const std::vector<SourcedData<DyeSeq, SourceCount<int>>>& dye_seqs,
        int num_timesteps,
        int num_channels,
        int radiometries_per_peptide,
        std::default_random_engine* generator,
        std::vector<SourcedData<Radiometry, SourceCount<int>>>* radiometries);

}  // namespace fluoroseq

#endif  // FLUOROSEQ_SIMULATION_GENERATE_RADIOMETRIES_H
