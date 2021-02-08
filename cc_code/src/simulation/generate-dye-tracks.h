/******************************************************************************\
* Author: Matthew Beauregard Smith                                             *
* Affiliation: The University of Texas at Austin                               *
* Department: Oden Institute and Institute for Cellular and Molecular Biology  *
* PI: Edward Marcotte                                                          *
* Project: Protein Fluorosequencing                                            *
\******************************************************************************/

#ifndef WHATPROT_SIMULATION_GENERATE_DYE_TRACKS_H
#define WHATPROT_SIMULATION_GENERATE_DYE_TRACKS_H

// Standard C++ library headers:
#include <random>
#include <unordered_map>
#include <vector>

// Local project headers:
#include "common/dye-seq.h"
#include "common/dye-track.h"
#include "common/error-model.h"
#include "common/sourced-data.h"

namespace whatprot {

void generate_dye_tracks(
        const ErrorModel& error_model,
        const std::vector<SourcedData<DyeSeq, SourceCount<int>>>& dye_seqs,
        int num_timesteps,
        int num_channels,
        int dye_tracks_per_peptide,
        std::default_random_engine* generator,
        std::vector<SourcedData<DyeTrack, SourceCount<int>>>* dye_tracks);

}  // namespace whatprot

#endif  // WHATPROT_SIMULATION_GENERATE_DYE_TRACKS_H
