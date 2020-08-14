// Author: Matthew Beauregard Smith (UT Austin)
//
// For MPI version, define compiler macro USE_MPI when building.
#include "dye_tracks_io.h"

#include <algorithm>
#include <fstream>
#include <string>
#include <vector>

#ifdef USE_MPI
#include <mpi.h>
#endif  // USE_MPI

#include "common/dye_track.h"
#include "common/sourced_data.h"

namespace fluoroseq {

namespace {
using std::copy;
using std::ifstream;
using std::string;
using std::vector;
}  // namespace

void read_dye_tracks(
        const string& filename,
        int* num_timesteps,
        int* num_channels,
        int* num_dye_tracks,
        SourcedData<DyeTrack*, SourceCountHitsList<int>*>*** dye_tracks) {
    int f_ints_size;
    int* f_ints;
    read_dye_tracks_raw(filename,
                        num_timesteps,
                        num_channels,
                        num_dye_tracks,
                        &f_ints_size,
                        &f_ints);
#ifdef USE_MPI
    broadcast_dye_tracks(num_timesteps,
                         num_channels,
                         num_dye_tracks,
                         &f_ints_size,
                         &f_ints);
#endif  // USE_MPI
    convert_dye_tracks_from_raw(*num_timesteps,
                                *num_channels,
                                *num_dye_tracks,
                                f_ints,
                                dye_tracks);
    delete[] f_ints;
}

void read_dye_tracks_raw(const string& filename,
                         int* num_timesteps,
                         int* num_channels,
                         int* num_dye_tracks,
                         int* f_ints_size,
                         int** f_ints) {
#ifdef USE_MPI
    int mpi_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
    if (mpi_rank != 0) {
        return;
    }
#endif  // USE_MPI
    ifstream f(filename);
    f >> *num_timesteps;
    f >> *num_channels;
    f >> *num_dye_tracks;
    vector<int> f_ints_vec;
    f_ints_vec.reserve((*num_dye_tracks)
                       * (3 + (*num_timesteps) * (*num_channels)));
    while (f.peek() != EOF) {
        int f_int;
        f >> f_int;
        f_ints_vec.push_back(f_int);
    }
    f.close();
    *f_ints_size = f_ints_vec.size();
    *f_ints = new int[*f_ints_size];
    copy(f_ints_vec.begin(), f_ints_vec.end(), *f_ints);
}

#ifdef USE_MPI
void broadcast_dye_tracks(int* num_timesteps,
                          int* num_channels,
                          int* num_dye_tracks,
                          int* f_ints_size,
                          int** f_ints) {
    int mpi_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
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
    MPI_Bcast(f_ints_size,
              1,  // count
              MPI_INT,
              0,  // root
              MPI_COMM_WORLD);
    if (mpi_rank != 0) {
        *f_ints = new int[*f_ints_size];
    }
    MPI_Bcast(*f_ints,
              *f_ints_size,
              MPI_INT,
              0,  // root
              MPI_COMM_WORLD);
}
#endif  // USE_MPI

void convert_dye_tracks_from_raw(
        int num_timesteps,
        int num_channels,
        int num_dye_tracks,
        int* f_ints,
        SourcedData<DyeTrack*, SourceCountHitsList<int>*>*** dye_tracks) {
    int index = 0;
    *dye_tracks = new SourcedData<
            DyeTrack*, SourceCountHitsList<int>*>*[num_dye_tracks];
    for (int i = 0; i < num_dye_tracks; i++) {
        DyeTrack* dye_track = new DyeTrack(num_timesteps, num_channels);
        for (int j = 0; j < num_timesteps * num_channels; j++) {
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
