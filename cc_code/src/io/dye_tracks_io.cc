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

void read_dye_tracks(
        const string& filename,
        int* num_timesteps,
        int* num_channels,
        int* num_dye_tracks,
        SourcedData<DyeTrack*,
                    SourceCountHitsList<int>*>*** dye_tracks) {
    ifstream f(filename);
    f >> *num_timesteps;
    f >> *num_channels;
    f >> *num_dye_tracks;
    *dye_tracks = new SourcedData<
            DyeTrack*, SourceCountHitsList<int>*>*[*num_dye_tracks];
    for (int i = 0; i < *num_dye_tracks; i++) {
        DyeTrack* dye_track = new DyeTrack(*num_timesteps, *num_channels);
        for (int j = 0; j < (*num_timesteps) * (*num_channels); j++) {
            f >> dye_track->counts[j];
        }
        int num_sources;
        f >> num_sources;
        SourceCountHits<int>** sources = new SourceCountHits<int>*[num_sources];
        for (int j = 0; j < num_sources; j++) {
            int id;
            f >> id;
            int count;
            f >> count;
            int hits;
            f >> hits;
            sources[j] = new SourceCountHits<int>(id, count, hits);
        }
        (*dye_tracks)[i] = new SourcedData<
                DyeTrack*, SourceCountHitsList<int>*>(
                        dye_track, new SourceCountHitsList<int>(num_sources,
                                                                sources));
    }
    f.close();
}

}  // namespace fluoroseq
