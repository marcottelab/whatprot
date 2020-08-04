// Author: Matthew Beauregard Smith (UT Austin)
#ifndef FLUOROSEQ_MPI_IO_DYE_TRACKS_MPI_IO_H
#define FLUOROSEQ_MPI_IO_DYE_TRACKS_MPI_IO_H

#include <string>

#include "common/dye_track.h"
#include "common/sourced_data.h"

namespace fluoroseq {

void mpi_read_dye_tracks(const std::string& filename,
                         int* num_timesteps,
                         int* num_channels,
                         int* num_dye_tracks,
                         SourcedData<DyeTrack*,
                                     SourceCountHitsList<int>*>*** dye_tracks);

void mpi_read_dye_tracks_master(
        const std::string& filename,
        int* num_timesteps,
        int* num_channels,
        int* num_dye_tracks,
        SourcedData<DyeTrack*, SourceCountHitsList<int>*>*** dye_tracks);

void mpi_read_dye_tracks_slave(
        const std::string& filename,
        int* num_timesteps,
        int* num_channels,
        int* num_dye_tracks,
        SourcedData<DyeTrack*, SourceCountHitsList<int>*>*** dye_tracks);

}  // namespace fluoroseq

#endif  // FLUOROSEQ_MPI_IO_DYE_TRACKS_MPI_IO_H