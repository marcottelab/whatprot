/******************************************************************************\
* Author: Matthew Beauregard Smith                                             *
* Affiliation: The University of Texas at Austin                               *
* Department: Oden Institute and Institute for Cellular and Molecular Biology  *
* PI: Edward Marcotte                                                          *
* Project: Protein Fluorosequencing                                            *
\******************************************************************************/

#ifndef WHATPROT_SIMULATION_GENERATE_RADIOMETRIES_H
#define WHATPROT_SIMULATION_GENERATE_RADIOMETRIES_H

// Standard C++ library headers:
#include <random>
#include <vector>

// Local project headers:
#include "common/dye-seq.h"
#include "common/radiometry.h"
#include "common/sourced-data.h"
#include "parameterization/model/sequencing-model.h"

namespace whatprot {

void generate_radiometries(
        const SequencingModel& seq_model,
        const std::vector<SourcedData<DyeSeq, SourceCount<int>>>& dye_seqs,
        int num_timesteps,
        int num_channels,
        int radiometries_per_peptide,
        std::default_random_engine* generator,
        std::vector<SourcedData<Radiometry, SourceCount<int>>>* radiometries);

}  // namespace whatprot

#endif  // WHATPROT_SIMULATION_GENERATE_RADIOMETRIES_H
