// Author: Matthew Beauregard Smith (UT Austin)
#ifndef FLUOROSEQ_IO_DYE_TRACKS_IO_H
#define FLUOROSEQ_IO_DYE_TRACKS_IO_H

#include <string>

#include "common/dye_track.h"
#include "common/sourced_data.h"

namespace fluoroseq {

void read_dye_tracks(const std::string& filename,
                     int* num_timesteps,
                     int* num_channels,
                     int* num_dye_tracks,
                     SourcedData<DyeTrack*,
                                 SourceCountMap<int>*>*** dye_tracks);

}  // namespace fluoroseq

#endif  // FLUOROSEQ_IO_DYE_TRACKS_IO_H