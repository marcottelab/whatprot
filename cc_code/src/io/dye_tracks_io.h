// Author: Matthew Beauregard Smith (UT Austin)
//
// For MPI version, define compiler macro USE_MPI when building.
#ifndef FLUOROSEQ_IO_DYE_TRACKS_IO_H
#define FLUOROSEQ_IO_DYE_TRACKS_IO_H

#include <string>
#include <vector>

#include "common/dye_track.h"
#include "common/sourced_data.h"

namespace fluoroseq {

void read_dye_tracks(
        const std::string& filename,
        int* num_timesteps,
        int* num_channels,
        std::vector<SourcedData<DyeTrack,
                                SourceCountHitsList<int>>>* dye_tracks);

void read_dye_tracks_raw(const std::string& filename,
                         int* num_timesteps,
                         int* num_channels,
                         int* num_dye_tracks,
                         int* f_ints_size,
                         int** f_ints);

#ifdef USE_MPI
void broadcast_dye_tracks(int* num_timesteps,
                          int* num_channels,
                          int* num_dye_tracks,
                          int* f_ints_size,
                          int** f_ints);
#endif  // USE_MPI

void convert_dye_tracks_from_raw(
        int num_timesteps,
        int num_channels,
        int num_dye_tracks,
        int* f_ints,
        std::vector<SourcedData<DyeTrack,
                                SourceCountHitsList<int>>>* dye_tracks);

void write_dye_tracks(
        const string& filename,
        int num_timesteps,
        int num_channels,
        const std::vector<SourcedData<DyeTrack,
                                      SourceCountHitsList<int>>>& dye_tracks);

// #ifdef USE_MPI

// void convert_raw_from_dye_tracks(
//         int num_timesteps,
//         int num_channels,
//         const std::vector<SourcedData<DyeTrack,
//                                       SourceCountHitsList<int>>>& dye_tracks,
//         int* f_ints_size,
//         int** f_ints);

// void gather_dye_tracks(int* num_dye_tracks, int* f_ints_size, int** f_ints);

// void merge_duplicate_dye_tracks(
//         std::vector<SourcedData<DyeTrack,
//                                 SourceCountHitsList<int>>>* dye_tracks);

// #endif  // USE_MPI

void write_dye_tracks_helper(
        const string& filename,
        int num_timesteps,
        int num_channels,
        const std::vector<SourcedData<DyeTrack,
                                      SourceCountHitsList<int>>>& dye_tracks);

}  // namespace fluoroseq

#endif  // FLUOROSEQ_IO_DYE_TRACKS_IO_H