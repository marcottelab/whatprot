// Author: Matthew Beauregard Smith (UT Austin)
#ifndef FLUOROSEQ_SIMULATION_GENERATE_RADIOMETRY_H
#define FLUOROSEQ_SIMULATION_GENERATE_RADIOMETRY_H

#include <random>

#include "common/dye_seq.h"
#include "common/radiometry.h"
#include "common/error_model.h"

namespace fluoroseq {

void generate_radiometry(const ErrorModel& error_model,
                         const DyeSeq& dye_seq,
                         int num_timesteps,
                         int num_channels,
                         std::default_random_engine* generator,
                         Radiometry* radiometry);

}  // namespace fluoroseq

#endif  // FLUOROSEQ_SIMULATION_GENERATE_RADIOMETRY_H
