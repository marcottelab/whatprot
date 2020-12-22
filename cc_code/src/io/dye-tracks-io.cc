/******************************************************************************\
* Author: Matthew Beauregard Smith                                             *
* Affiliation: The University of Texas at Austin                               *
* Department: Oden Institute and Institute for Cellular and Molecular Biology  *
* PI: Edward Marcotte                                                          *
* Project: Protein Fluorosequencing                                            *
\******************************************************************************/

// Defining symbols from header:
#include "dye-tracks-io.h"

// Standard C++ library headers:
#include <algorithm>
#include <fstream>
#include <string>
#include <vector>

// Local project headers:
#include "common/dye-track.h"
#include "common/sourced-data.h"

namespace fluoroseq {

namespace {
using std::copy;
using std::ifstream;
using std::ofstream;
using std::string;
using std::vector;
}  // namespace

void read_dye_tracks(
        const string& filename,
        int* num_timesteps,
        int* num_channels,
        vector<SourcedData<DyeTrack, SourceCountHitsList<int>>>* dye_tracks) {
    int num_dye_tracks;
    int f_ints_size;
    int* f_ints;
    read_dye_tracks_raw(filename,
                        num_timesteps,
                        num_channels,
                        &num_dye_tracks,
                        &f_ints_size,
                        &f_ints);
    convert_dye_tracks_from_raw(
            *num_timesteps, *num_channels, num_dye_tracks, f_ints, dye_tracks);
    delete[] f_ints;
}

void read_dye_tracks_raw(const string& filename,
                         int* num_timesteps,
                         int* num_channels,
                         int* num_dye_tracks,
                         int* f_ints_size,
                         int** f_ints) {
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

void convert_dye_tracks_from_raw(
        int num_timesteps,
        int num_channels,
        int num_dye_tracks,
        int* f_ints,
        vector<SourcedData<DyeTrack, SourceCountHitsList<int>>>* dye_tracks) {
    int index = 0;
    dye_tracks->reserve(num_dye_tracks);
    for (int i = 0; i < num_dye_tracks; i++) {
        DyeTrack dye_track(num_timesteps, num_channels);
        for (int j = 0; j < num_timesteps * num_channels; j++) {
            dye_track.counts[j] = f_ints[index];
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
        dye_tracks->push_back(SourcedData<DyeTrack, SourceCountHitsList<int>>(
                dye_track, SourceCountHitsList<int>(num_sources, sources)));
    }
}

void write_dye_tracks(
        const string& filename,
        int num_timesteps,
        int num_channels,
        const vector<SourcedData<DyeTrack, SourceCountHitsList<int>>>&
                dye_tracks) {
    const vector<SourcedData<DyeTrack, SourceCountHitsList<int>>>&
            new_dye_tracks = dye_tracks;
    write_dye_tracks_helper(
            filename, num_timesteps, num_channels, new_dye_tracks);
}

void write_dye_tracks_helper(
        const string& filename,
        int num_timesteps,
        int num_channels,
        const vector<SourcedData<DyeTrack, SourceCountHitsList<int>>>&
                dye_tracks) {
    ofstream f(filename);
    f << num_timesteps << "\n";
    f << num_channels << "\n";
    f << dye_tracks.size() << "\n";
    for (const auto& dye_track : dye_tracks) {
        for (int i = 0; i < num_timesteps * num_channels; i++) {
            f << dye_track.value.counts[i] << "\t";
        }
        f << dye_track.source.num_sources << "\t";
        for (int i = 0; i < dye_track.source.num_sources; i++) {
            f << dye_track.source.sources[i]->source << "\t";
            f << dye_track.source.sources[i]->count << "\t";
            f << dye_track.source.sources[i]->hits;
            if (i < dye_track.source.num_sources - 1) {
                f << "\t";
            }
        }
        f << "\n";
    }
    f.close();
}

}  // namespace fluoroseq
