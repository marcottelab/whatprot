// Author: Matthew Beauregard Smith (UT Austin)
#ifndef FLUOROSEQ_SIMULATION_GENERATE_DYE_TRACK_H
#define FLUOROSEQ_SIMULATION_GENERATE_DYE_TRACK_H

#include <random>

#include "common/dye_seq.h"
#include "common/dye_track.h"
#include "common/error_model.h"

namespace fluoroseq {

void generate_dye_track(const ErrorModel& error_model,
                        const DyeSeq& dye_seq,
                        int num_timesteps,
                        int num_channels,
                        std::default_random_engine* generator,
                        DyeTrack* dye_track);

}  // namespace fluoroseq

#endif  // FLUOROSEQ_SIMULATION_GENERATE_DYE_TRACK_H
