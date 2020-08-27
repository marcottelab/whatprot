// Author: Matthew Beauregard Smith (UT Austin)
#include "generate_radiometry.h"

#include <cmath>
#include <random>

#include "simulation/generate_dye_track.h"

namespace fluoroseq {

namespace {
using std::default_random_engine;
using std::log;
using std::lognormal_distribution;
}  // namespace

void generate_radiometry(const ErrorModel& error_model,
                         const DyeSeq& dye_seq,
                         int num_timesteps,
                         int num_channels,
                         default_random_engine* generator,
                         Radiometry* radiometry) {
    DyeTrack dye_track(num_timesteps, num_channels);
    generate_dye_track(error_model,
                       dye_seq,
                       num_timesteps,
                       num_channels,
                       generator,
                       &dye_track);
    for (int t = 0; t < num_timesteps; t++) {
        for (int c = 0; c < num_channels; c++) {
            lognormal_distribution<double> lognormal(
                    log(error_model.mu * (double) dye_track(t, c)),
                    error_model.sigma);
            (*radiometry)(t, c) = lognormal(*generator);
        }
    }
}

}  // namespace fluoroseq
