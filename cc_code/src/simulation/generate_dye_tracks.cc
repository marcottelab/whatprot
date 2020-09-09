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
        vector<SourcedData<DyeTrack, SourceCountHitsList<int>>>* dye_tracks) {
    // Maps from DyeTrack to a map from the source to the number of hits.
    unordered_map<DyeTrack, unordered_map<int, int>> dye_track_map;
    dye_track_map.reserve(dye_seqs.size() * dye_tracks_per_dye_seq);
    unordered_map<int, int> source_count_map;
    source_count_map.reserve(dye_seqs.size());
    for (const SourcedData<DyeSeq, SourceCount<int>>& dye_seq : dye_seqs) {
        int source = dye_seq.source.source;
        int count = dye_seq.source.count;
        source_count_map[source] = count;
        for (int j = 0; j < dye_tracks_per_dye_seq; j++) {
            DyeTrack dye_track(num_timesteps, num_channels);
            generate_dye_track(error_model,
                               dye_seq.value,
                               num_timesteps,
                               num_channels,
                               generator,
                               &dye_track);
            unordered_map<int, int>& source_hit_map = dye_track_map[
                    move(dye_track)];
            // Set to 0 if source is not in this source_hit_map
            // Important question: is this necessary, or does the value default
            // to 0?
            if (source_hit_map.find(source) == source_hit_map.end()) {
                source_hit_map[source] = 0;
            }
            source_hit_map[source]++;
        }
    }
    dye_tracks->reserve(dye_track_map.size());
    for (const auto& dye_track_entry : dye_track_map) {
        const DyeTrack& dye_track = dye_track_entry.first;
        const unordered_map<int, int>& source_hits_map = dye_track_entry.second;
        int num_sources = source_hits_map.size();
        SourceCountHits<int>** sources = new SourceCountHits<int>*[num_sources];
        int index = 0;
        for (const auto& source_hits : source_hits_map) {
            int id = source_hits.first;
            int count = source_count_map[id];
            int hits = source_hits.second;
            sources[index] = new SourceCountHits<int>(id, count, hits);
            index++;
        }
        dye_tracks->push_back(SourcedData<DyeTrack, SourceCountHitsList<int>>(
                dye_track,
                move(SourceCountHitsList<int>(num_sources, sources))));
    }
}

}  // namespace fluoroseq
