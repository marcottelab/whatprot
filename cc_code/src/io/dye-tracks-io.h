/******************************************************************************\
* Author: Matthew Beauregard Smith                                             *
* Affiliation: The University of Texas at Austin                               *
* Department: Oden Institute and Institute for Cellular and Molecular Biology  *
* PI: Edward Marcotte                                                          *
* Project: Protein Fluorosequencing                                            *
\******************************************************************************/

#ifndef FLUOROSEQ_IO_DYE_TRACKS_IO_H
#define FLUOROSEQ_IO_DYE_TRACKS_IO_H

// Standard C++ library headers:
#include <string>
#include <vector>

// Local project headers:
#include "common/dye-track.h"
#include "common/sourced-data.h"

namespace fluoroseq {

void read_dye_tracks(
        const std::string& filename,
        int* num_timesteps,
        int* num_channels,
        std::vector<SourcedData<DyeTrack, SourceCountHitsList<int>>>*
                dye_tracks);

void read_dye_tracks_raw(const std::string& filename,
                         int* num_timesteps,
                         int* num_channels,
                         int* num_dye_tracks,
                         int* f_ints_size,
                         int** f_ints);

void convert_dye_tracks_from_raw(
        int num_timesteps,
        int num_channels,
        int num_dye_tracks,
        int* f_ints,
        std::vector<SourcedData<DyeTrack, SourceCountHitsList<int>>>*
                dye_tracks);

void write_dye_tracks(
        const string& filename,
        int num_timesteps,
        int num_channels,
        const std::vector<SourcedData<DyeTrack, SourceCountHitsList<int>>>&
                dye_tracks);

void write_dye_tracks_helper(
        const string& filename,
        int num_timesteps,
        int num_channels,
        const std::vector<SourcedData<DyeTrack, SourceCountHitsList<int>>>&
                dye_tracks);

}  // namespace fluoroseq

#endif  // FLUOROSEQ_IO_DYE_TRACKS_IO_H