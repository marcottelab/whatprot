// Author: Matthew Beauregard Smith (UT Austin)
#ifndef FLUOROSEQ_SIMULATION_GENERATE_DYE_TRACKS_H
#define FLUOROSEQ_SIMULATION_GENERATE_DYE_TRACKS_H

#include <random>
#include <unordered_map>
#include <vector>

#include "common/dye_seq.h"
#include "common/dye_track.h"
#include "common/error_model.h"
#include "common/sourced_data.h"

namespace fluoroseq {

void generate_dye_tracks(
        const ErrorModel& error_model,
        const std::vector<SourcedData<DyeSeq, SourceCount<int>>>& dye_seqs,
        int num_timesteps,
        int num_channels,
        int dye_tracks_per_dye_seq,
        std::default_random_engine* generator,
        std::vector<SourcedData<DyeTrack, SourceCount<int>>>* dye_tracks);

}  // namespace fluoroseq

#endif  // FLUOROSEQ_SIMULATION_GENERATE_DYE_TRACKS_H
