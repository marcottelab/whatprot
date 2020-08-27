// Author: Matthew Beauregard Smith (UT Austin)
#ifndef FLUOROSEQ_SIMULATION_GENERATE_RADIOMETRIES_H
#define FLUOROSEQ_SIMULATION_GENERATE_RADIOMETRIES_H

#include <random>
#include <vector>

#include "common/dye_seq.h"
#include "common/radiometry.h"
#include "common/error_model.h"
#include "common/sourced_data.h"

namespace fluoroseq {

void generate_radiometries(
        const ErrorModel& error_model,
        const std::vector<SourcedData<DyeSeq, SourceCount<int>>>& dye_seqs,
        int num_timesteps,
        int num_channels,
        int radiometries_per_dye_seq,
        std::default_random_engine* generator,
        std::vector<SourcedData<Radiometry, SourceCount<int>>>* radiometries);

}  // namespace fluoroseq

#endif  // FLUOROSEQ_SIMULATION_GENERATE_RADIOMETRIES_H
