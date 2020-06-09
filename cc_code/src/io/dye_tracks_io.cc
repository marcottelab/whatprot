// Author: Matthew Beauregard Smith (UT Austin)
#include "dye_tracks_io.h"

#include <fstream>
#include <string>

#include "common/dye_track.h"

namespace fluoroseq {

namespace {
using std::ifstream;
using std::string;
}  // namespace

}  // namespace fluoroseq
#include <iostream>
namespace fluoroseq {

void read_dye_tracks(const string& filename,
                     int* num_timesteps,
                     int* num_channels,
                     int* num_dye_tracks,
                     SourcedData<DyeTrack*,
                                 SourceCountMap<int>*>*** dye_tracks) {
    ifstream f(filename);
    f >> *num_timesteps;
    f >> *num_channels;
    f >> *num_dye_tracks;
    // std::cout << "CHECKPOINT ONE\n";
    *dye_tracks = new SourcedData<DyeTrack*,
                                  SourceCountMap<int>*>*[*num_dye_tracks];
    // std::cout << "CHECKPOINT TWO\n";
    for (int i = 0; i < *num_dye_tracks; i++) {
        // std::cout << "iteration " << i << "\n";
        // std::cout << "SPOT ALPHA\n";
        DyeTrack* dye_track = new DyeTrack(*num_timesteps, *num_channels);
        // std::cout << "SPOT BETA\n";
        for (int j = 0; j < (*num_timesteps) * (*num_channels); j++) {
            f >> dye_track->counts[j];
        }
        // std::cout << "SPOT DELTA\n";
        int num_sources;
        f >> num_sources;
        SourceWithCount<int>** sources = new SourceWithCount<int>*[num_sources];
        // std::cout << "SPOT EPSILON\n";
        for (int j = 0; j < num_sources; j++) {
            // std::cout << "j is " << j << "\n";
            int id;
            f >> id;
            // std::cout << "id is " << id << "\n";
            int count;
            f >> count;
            // std::cout << "count is " << count << "\n";
            // std::cout << "a\n";
            sources[j] = new SourceWithCount<int>(id, count);
            // std::cout << "b\n";
        }
        // std::cout << "SPOT GAMMA\n";
        (*dye_tracks)[i] = new SourcedData<DyeTrack*, SourceCountMap<int>*>(
                dye_track, new SourceCountMap<int>(num_sources, sources));
    }
    // std::cout << "CHECKPOINT THREE\n";
    f.close();
    // std::cout << "CHECKPOINT FOUR\n";
}

}  // namespace fluoroseq
