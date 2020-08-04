// Author: Matthew Beauregard Smith (UT Austin)
#include "dye_tracks_mpi_io.h"

#include <fstream>
#include <string>
#include <vector>

#include <mpi.h>

#include "common/dye_track.h"
#include "common/sourced_data.h"

namespace fluoroseq {

namespace {
using std::ifstream;
using std::string;
using std::vector;
}  // namespace

void mpi_read_dye_tracks(const std::string& filename,
                         int* num_timesteps,
                         int* num_channels,
                         int* num_dye_tracks,
                         SourcedData<DyeTrack*,
                                     SourceCountHitsList<int>*>*** dye_tracks) {
    int mpi_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
    if (mpi_rank == 0) {
        mpi_read_dye_tracks_master(filename,
                                   num_timesteps,
                                   num_channels,
                                   num_dye_tracks,
                                   dye_tracks);
    } else {
        mpi_read_dye_tracks_slave(filename,
                                  num_timesteps,
                                  num_channels,
                                  num_dye_tracks,
                                  dye_tracks);
    }
}

void mpi_read_dye_tracks_master(
        const std::string& filename,
        int* num_timesteps,
        int* num_channels,
        int* num_dye_tracks,
        SourcedData<DyeTrack*, SourceCountHitsList<int>*>*** dye_tracks) {
    ifstream f(filename);
    f >> *num_timesteps;
    MPI_Bcast(num_timesteps,
              1,  // count
              MPI_INT,
              0,  // root
              MPI_COMM_WORLD);
    f >> *num_channels;
    MPI_Bcast(num_channels,
              1,  // count
              MPI_INT,
              0,  // root
              MPI_COMM_WORLD);
    f >> *num_dye_tracks;
    MPI_Bcast(num_dye_tracks,
              1,  // count
              MPI_INT,
              0,  // root
              MPI_COMM_WORLD);
    vector<int> f_ints;
    f_ints.reserve((*num_dye_tracks)
                   * (3 + (*num_timesteps) * (*num_channels)));
    while (f.peek() != EOF) {
        int f_int;
        f >> f_int;
        f_ints.push_back(f_int);
    }
    f.close();
    int f_ints_size = f_ints.size();
    MPI_Bcast(&f_ints_size,
              1,  // count
              MPI_INT,
              0,  // root
              MPI_COMM_WORLD);
    MPI_Bcast(&f_ints.front(),
              f_ints.size(),
              MPI_INT,
              0,  // root
              MPI_COMM_WORLD);
    int index = 0;
    *dye_tracks = new SourcedData<
            DyeTrack*, SourceCountHitsList<int>*>*[*num_dye_tracks];
    for (int i = 0; i < *num_dye_tracks; i++) {
        DyeTrack* dye_track = new DyeTrack(*num_timesteps, *num_channels);
        for (int j = 0; j < (*num_timesteps) * (*num_channels); j++) {

            dye_track->counts[j] = f_ints[index];
            index++;
        }
        int num_sources = f_ints[index];
        index++;
        SourceCountHits<int>** sources = new SourceCountHits<int>*[num_sources];
        for (int j = 0; j < num_sources; j++) {
            int id = f_ints[index];
            index++;
            int count = f_ints[index];
            index++;
            int hits = f_ints[index];
            index++;
            sources[j] = new SourceCountHits<int>(id, count, hits);
        }
        (*dye_tracks)[i] = new SourcedData<
                DyeTrack*, SourceCountHitsList<int>*>(
                        dye_track,
                        new SourceCountHitsList<int>(num_sources, sources));
    }
}

void mpi_read_dye_tracks_slave(
        const std::string& filename,
        int* num_timesteps,
        int* num_channels,
        int* num_dye_tracks,
        SourcedData<DyeTrack*, SourceCountHitsList<int>*>*** dye_tracks) {
    MPI_Bcast(num_timesteps,
              1,  // count
              MPI_INT,
              0,  // root
              MPI_COMM_WORLD);
    MPI_Bcast(num_channels,
              1,  // count
              MPI_INT,
              0,  // root
              MPI_COMM_WORLD);
    MPI_Bcast(num_dye_tracks,
              1,  // count
              MPI_INT,
              0,  // root
              MPI_COMM_WORLD);
    int f_ints_size;
    MPI_Bcast(&f_ints_size,
              1,  // count
              MPI_INT,
              0,  // root
              MPI_COMM_WORLD);
    int* f_ints = new int[f_ints_size];
    MPI_Bcast(f_ints,
              f_ints_size,
              MPI_INT,
              0,  // root
              MPI_COMM_WORLD);
    int index = 0;
    *dye_tracks = new SourcedData<
            DyeTrack*, SourceCountHitsList<int>*>*[*num_dye_tracks];
    for (int i = 0; i < *num_dye_tracks; i++) {
        DyeTrack* dye_track = new DyeTrack(*num_timesteps, *num_channels);
        for (int j = 0; j < (*num_timesteps) * (*num_channels); j++) {

            dye_track->counts[j] = f_ints[index];
            index++;
        }
        int num_sources = f_ints[index];
        index++;
        SourceCountHits<int>** sources = new SourceCountHits<int>*[num_sources];
        for (int j = 0; j < num_sources; j++) {
            int id = f_ints[index];
            index++;
            int count = f_ints[index];
            index++;
            int hits = f_ints[index];
            index++;
            sources[j] = new SourceCountHits<int>(id, count, hits);
        }
        (*dye_tracks)[i] = new SourcedData<
                DyeTrack*, SourceCountHitsList<int>*>(
                        dye_track,
                        new SourceCountHitsList<int>(num_sources, sources));
    }
}

}  // namespace fluoroseq
