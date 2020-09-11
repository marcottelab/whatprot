// Author: Matthew Beauregard Smith (UT Austin)
//
// For MPI version, define compiler macro USE_MPI when building.
#include "dedup_dye_tracks.h"

#ifdef USE_MPI

namespace fluoroseq {

void dedup_dye_tracks(
        vector<SourcedData<DyeTrack, SourceCount<int>>>* dye_tracks_in,
        vector<SourcedData<DyeTrack,
                           SourceCountHitsList<int>>>* dye_tracks_out) {}

}  // namespace fluoroseq

#else  // USE_MPI

#include <unordered_map>
#include <utility>  // for std::move
#include <vector>

#include "common/dye_track.h"
#include "common/sourced_data.h"

namespace fluoroseq {

namespace {
using std::move;
using std::unordered_map;
using std::vector;
}  // namespace

void dedup_dye_tracks(
        vector<SourcedData<DyeTrack, SourceCount<int>>>* dye_tracks_in,
        vector<SourcedData<DyeTrack,
                           SourceCountHitsList<int>>>* dye_tracks_out) {
    unordered_map<DyeTrack,
                  unordered_map<int, SourceCountHits<int>>> dye_track_map;
    dye_track_map.reserve(dye_tracks_in->size());
    while (!dye_tracks_in->empty()) {
        SourcedData<DyeTrack,
                    SourceCount<int>>& dye_track = dye_tracks_in->back();
        int source = dye_track.source.source;
        int count = dye_track.source.count;
        unordered_map<int, SourceCountHits<int>>& source_map = dye_track_map[
                move(dye_track.value)];
        if (source_map.find(source) == source_map.end()) {
            source_map[source] = SourceCountHits<int>(source, count, 1);
        } else {
            source_map[source].hits++;
        }
        dye_tracks_in->pop_back();
    }
    dye_tracks_out->reserve(dye_track_map.size());
    for (const auto& dye_track_map_entry : dye_track_map) {
        const DyeTrack& dye_track = dye_track_map_entry.first;
        const unordered_map<int, SourceCountHits<
                int>>& source_map = dye_track_map_entry.second;
        int num_sources = source_map.size();
        SourceCountHits<int>** sources = new SourceCountHits<int>*[num_sources];
        int index = 0;
        for (const auto& source_map_entry : source_map) {
            const SourceCountHits<
                    int>& source_count_hits = source_map_entry.second;
            sources[index] = new SourceCountHits<int>(source_count_hits);
            index++;
        }
        dye_tracks_out->push_back(
                move(SourcedData<DyeTrack, SourceCountHitsList<int>>(
                        dye_track,
                        move(SourceCountHitsList<int>(num_sources, sources)))));
    }
}

}  // namespace fluoroseq

#endif  // USE_MPI
