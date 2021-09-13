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
using std::normal_distribution;
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
                int num_dyes = dye_track(t, c);
                double mu = seq_model.channel_models[c]->mu * num_dyes;
                double sig = seq_model.channel_models[c]->sigma(num_dyes);
                normal_distribution<double> normal(mu, sig);
                (*radiometry)(t, c) = normal(*generator);
        }
    }
}

}  // namespace whatprot
