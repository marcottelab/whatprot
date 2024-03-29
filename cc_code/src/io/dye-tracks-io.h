/******************************************************************************\
* Author: Matthew Beauregard Smith                                             *
* Affiliation: The University of Texas at Austin                               *
* Department: Oden Institute and Institute for Cellular and Molecular Biology  *
* PI: Edward Marcotte                                                          *
* Project: Protein Fluorosequencing                                            *
\******************************************************************************/

#ifndef WHATPROT_IO_DYE_TRACKS_IO_H
#define WHATPROT_IO_DYE_TRACKS_IO_H

// Standard C++ library headers:
#include <string>
#include <vector>

// Local project headers:
#include "common/dye-track.h"
#include "common/sourced-data.h"

namespace whatprot {

void read_dye_tracks(
        const std::string& filename,
        unsigned int* num_timesteps,
        unsigned int* num_channels,
        std::vector<SourcedData<DyeTrack, SourceCountHitsList<int>>>*
                dye_tracks);

void read_dye_tracks_raw(const std::string& filename,
                         unsigned int* num_timesteps,
                         unsigned int* num_channels,
                         unsigned int* num_dye_tracks,
                         unsigned int* f_ints_size,
                         unsigned int** f_ints);

void convert_dye_tracks_from_raw(
        unsigned int num_timesteps,
        unsigned int num_channels,
        unsigned int num_dye_tracks,
        unsigned int* f_ints,
        std::vector<SourcedData<DyeTrack, SourceCountHitsList<int>>>*
                dye_tracks);

void write_dye_tracks(
        const string& filename,
        unsigned int num_timesteps,
        unsigned int num_channels,
        const std::vector<SourcedData<DyeTrack, SourceCountHitsList<int>>>&
                dye_tracks);

void write_dye_tracks_helper(
        const string& filename,
        unsigned int num_timesteps,
        unsigned int num_channels,
        const std::vector<SourcedData<DyeTrack, SourceCountHitsList<int>>>&
                dye_tracks);

}  // namespace whatprot

#endif  // WHATPROT_IO_DYE_TRACKS_IO_H