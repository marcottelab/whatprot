// Author: Matthew Beauregard Smith (UT Austin)
//
// For MPI version, define compiler macro USE_MPI when building.
#ifndef FLUOROSEQ_SIMULATION_DEDUP_DYE_TRACKS_H
#define FLUOROSEQ_SIMULATION_DEDUP_DYE_TRACKS_H

#include <unordered_map>
#include <vector>

#include "common/dye_track.h"
#include "common/sourced_data.h"
#include "keyvalue.h"

namespace fluoroseq {

void dedup_dye_tracks(
        int num_timesteps,
        int num_channels,
        std::vector<SourcedData<DyeTrack, SourceCount<int>>>* dye_tracks_in,
        std::vector<SourcedData<DyeTrack,
                                SourceCountHitsList<int>>>* dye_tracks_out);

void map_dye_tracks(int index, MAPREDUCE_NS::KeyValue* kv, void* ptr);

void reduce_dye_tracks(char* key,
                       int keybytes,
                       char* multivalue,
                       int nvalues,
                       int* valuebytes,
                       MAPREDUCE_NS::KeyValue* kv,
                       void* ptr);

void reduce_dye_tracks_helper(
        char* multivalue,
        int nvalues,
        std::unordered_map<int, SourceCountHits<int>>* source_map);

class OutputInfo {
public:
    OutputInfo(int num_timesteps,
               int num_channels,
               std::vector<SourcedData<
                       DyeTrack, SourceCountHitsList<int>>>* dye_tracks_out);
    int num_timesteps;
    int num_channels;
    std::vector<
            SourcedData<DyeTrack,
                        SourceCountHitsList<int>>>* dye_tracks_out;  // not ownd
};

void output_dye_tracks(char* key,
                       int keybytes,
                       char* value,
                       int valuebytes,
                       void* ptr);

}  // namespace fluoroseq

#endif  // FLUOROSEQ_SIMULATION_DEDUP_DYE_TRACKS_H
