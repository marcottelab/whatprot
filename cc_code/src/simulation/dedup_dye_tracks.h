// Author: Matthew Beauregard Smith (UT Austin)
//
// For MPI version, define compiler macro USE_MPI when building.
#ifndef FLUOROSEQ_SIMULATION_DEDUP_DYE_TRACKS_H
#define FLUOROSEQ_SIMULATION_DEDUP_DYE_TRACKS_H

#include <vector>

#include "common/dye_track.h"
#include "common/sourced_data.h"

namespace fluoroseq {

void dedup_dye_tracks(
        std::vector<SourcedData<DyeTrack, SourceCount<int>>>* dye_tracks_in,
        std::vector<SourcedData<DyeTrack,
                                SourceCountHitsList<int>>>* dye_tracks_out);

}  // namespace fluoroseq

#endif  // FLUOROSEQ_SIMULATION_DEDUP_DYE_TRACKS_H
