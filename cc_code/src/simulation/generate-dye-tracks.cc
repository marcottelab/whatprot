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
using std::move;
using std::vector;
}  // namespace

void generate_dye_tracks(
        const SequencingModel& seq_model,
        const vector<SourcedData<DyeSeq, SourceCount<int>>>& dye_seqs,
        unsigned int num_timesteps,
        unsigned int num_channels,
        unsigned int dye_tracks_per_peptide,
        default_random_engine* generator,
        vector<SourcedData<DyeTrack, SourceCount<int>>>* dye_tracks) {
    dye_tracks->reserve(dye_seqs.size());
    for (const SourcedData<DyeSeq, SourceCount<int>>& dye_seq : dye_seqs) {
        // We want to generate a certain number of radiometries per peptide,
        // not per dye_seq. Therefore we do this on repeat for each peptide
        // that produced this dye_seq.
        for (int i = 0; i < dye_seq.source.count; i++) {
            for (unsigned int j = 0; j < dye_tracks_per_peptide; j++) {
                DyeTrack dye_track(num_timesteps, num_channels);
                generate_dye_track(seq_model,
                                   dye_seq.value,
                                   num_timesteps,
                                   num_channels,
                                   generator,
                                   &dye_track);
                // Ignore any DyeTrack with all 0s because it wouldn't be
                // detectable. Any DyeTrack with all 0s at the 0th timestep will
                // have all 0s throughout.
                bool nontrivial = false;
                for (unsigned int c = 0; c < num_channels; c++) {
                    if (dye_track(0, c) != 0) {
                        nontrivial = true;
                    }
                }
                if (nontrivial) {
                    dye_tracks->push_back(
                            move(SourcedData<DyeTrack, SourceCount<int>>(
                                    move(dye_track), dye_seq.source)));
                }
            }
        }
    }
}

}  // namespace whatprot
