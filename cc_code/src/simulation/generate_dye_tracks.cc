// Author: Matthew Beauregard Smith (UT Austin)
#include "generate_dye_tracks.h"

#include <random>
#include <unordered_map>
#include <utility>  // for std::move
#include <vector>

#include "common/dye_seq.h"
#include "common/dye_track.h"
#include "common/error_model.h"
#include "common/sourced_data.h"
#include "simulation/generate_dye_track.h"

namespace fluoroseq {

namespace {
using std::default_random_engine;
using std::move;
using std::unordered_map;
using std::vector;
}  // namespace

void generate_dye_tracks(
        const ErrorModel& error_model,
        const vector<SourcedData<DyeSeq, SourceCount<int>>>& dye_seqs,
        int num_timesteps,
        int num_channels,
        int dye_tracks_per_dye_seq,
        default_random_engine* generator,
        vector<SourcedData<DyeTrack, SourceCount<int>>>* dye_tracks) {
    dye_tracks->reserve(dye_seqs.size());
    for (const SourcedData<DyeSeq, SourceCount<int>>& dye_seq : dye_seqs) {
        for (int j = 0; j < dye_tracks_per_dye_seq; j++) {
            DyeTrack dye_track(num_timesteps, num_channels);
            generate_dye_track(error_model,
                               dye_seq.value,
                               num_timesteps,
                               num_channels,
                               generator,
                               &dye_track);
            dye_tracks->push_back(
                    move(SourcedData<DyeTrack, SourceCount<int>>(
                            move(dye_track), dye_seq.source)));
        }
    }
}

}  // namespace fluoroseq
