/******************************************************************************\
* Author: Matthew Beauregard Smith                                             *
* Affiliation: The University of Texas at Austin                               *
* Department: Oden Institute and Institute for Cellular and Molecular Biology  *
* PI: Edward Marcotte                                                          *
* Project: Protein Fluorosequencing                                            *
\******************************************************************************/

// Defining symbols from header:
#include "generate-dye-tracks.h"

// Standard C++ library headers:
#include <random>
#include <utility>  // for std::move
#include <vector>

// Local project headers:
#include "common/dye-seq.h"
#include "common/dye-track.h"
#include "common/sourced-data.h"
#include "parameterization/model/sequencing-model.h"
#include "simulation/generate-dye-track.h"

namespace whatprot {

namespace {
using std::default_random_engine;
using std::discrete_distribution;
using std::move;
using std::vector;
}  // namespace

void generate_dye_tracks(
        const SequencingModel& seq_model,
        const vector<SourcedData<DyeSeq, SourceCount<int>>>& dye_seqs,
        unsigned int num_timesteps,
        unsigned int num_channels,
        unsigned int num_to_generate,
        default_random_engine* generator,
        vector<SourcedData<DyeTrack, SourceCount<int>>>* dye_tracks) {
    // We want the dye tracks generated based on a uniform distribution of
    // peptides, not dye-seqs. We therefore need a discrete_distribution,
    // because it is weighted.
    vector<double> index_to_weight(dye_seqs.size(), 0);
    for (unsigned int i = 0; i < dye_seqs.size(); i++) {
        index_to_weight[i] = (double)dye_seqs[i].source.count;
    }
    discrete_distribution<unsigned int> random_dye_seq_idx(
            index_to_weight.begin(), index_to_weight.end());
    for (unsigned int i = 0; i < num_to_generate; i++) {
        unsigned int dye_seq_idx = random_dye_seq_idx(*generator);
        DyeTrack dye_track(num_timesteps, num_channels);
        generate_dye_track(seq_model,
                           dye_seqs[dye_seq_idx].value,
                           num_timesteps,
                           num_channels,
                           generator,
                           &dye_track);
        // Ignore any DyeTrack with all 0s because it wouldn't be detectable.
        // Any DyeTrack with all 0s at the 0th timestep will have all 0s
        // throughout.
        bool trivial = true;
        for (unsigned int c = 0; c < num_channels; c++) {
            if (dye_track(0, c) != 0) {
                trivial = false;
            }
        }
        if (!trivial) {
            dye_tracks->push_back(move(SourcedData<DyeTrack, SourceCount<int>>(
                    move(dye_track), dye_seqs[dye_seq_idx].source)));
        }
    }
}

}  // namespace whatprot
