/******************************************************************************\
* Author: Matthew Beauregard Smith                                             *
* Affiliation: The University of Texas at Austin                               *
* Department: Oden Institute and Institute for Cellular and Molecular Biology  *
* PI: Edward Marcotte                                                          *
* Project: Protein Fluorosequencing                                            *
\******************************************************************************/

// Defining symbols from header:
#include "generate-radiometry.h"

// Standard C++ library headers:
#include <cmath>
#include <random>

// Local project headers:
#include "parameterization/model/sequencing-model.h"
#include "simulation/generate-dye-track.h"

namespace whatprot {

namespace {
using std::default_random_engine;
using std::log;
using std::lognormal_distribution;
}  // namespace

void generate_radiometry(const SequencingModel& seq_model,
                         const DyeSeq& dye_seq,
                         unsigned int num_timesteps,
                         unsigned int num_channels,
                         default_random_engine* generator,
                         Radiometry* radiometry) {
    DyeTrack dye_track(num_timesteps, num_channels);
    generate_dye_track(seq_model,
                       dye_seq,
                       num_timesteps,
                       num_channels,
                       generator,
                       &dye_track);
    for (unsigned int t = 0; t < num_timesteps; t++) {
        for (unsigned int c = 0; c < num_channels; c++) {
            if (dye_track(t, c) > 0) {
                double mu = seq_model.channel_models[c]->mu;
                double sigma = seq_model.channel_models[c]->sigma;
                lognormal_distribution<double> lognormal(
                        mu + log((double)dye_track(t, c)), sigma);
                (*radiometry)(t, c) = lognormal(*generator);
            } else {
                (*radiometry)(t, c) = 0.0;
            }
        }
    }
}

}  // namespace whatprot
