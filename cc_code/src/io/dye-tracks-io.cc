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

namespace whatprot {

namespace {
using std::copy;
using std::ifstream;
using std::ofstream;
using std::string;
using std::vector;
}  // namespace

void read_dye_tracks(
        const string& filename,
        unsigned int* num_timesteps,
        unsigned int* num_channels,
        vector<SourcedData<DyeTrack, SourceCountHitsList<int>>>* dye_tracks) {
    unsigned int num_dye_tracks;
    ifstream f(filename);
    f >> *num_timesteps;
    f >> *num_channels;
    f >> num_dye_tracks;
    dye_tracks->reserve(num_dye_tracks);
    for (unsigned int i = 0; i < num_dye_tracks; i++) {
        DyeTrack dye_track(*num_timesteps, *num_channels);
        for (unsigned int j = 0; j < (*num_timesteps) * (*num_channels); j++) {
            f >> dye_track.counts[j];
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
        dye_tracks->push_back(SourcedData<DyeTrack, SourceCountHitsList<int>>(
                dye_track, SourceCountHitsList<int>(num_sources, sources)));
    }
    f.close();
}

void write_dye_tracks(
        const string& filename,
        unsigned int num_timesteps,
        unsigned int num_channels,
        const vector<SourcedData<DyeTrack, SourceCountHitsList<int>>>&
                dye_tracks) {
    const vector<SourcedData<DyeTrack, SourceCountHitsList<int>>>&
            new_dye_tracks = dye_tracks;
    write_dye_tracks_helper(
            filename, num_timesteps, num_channels, new_dye_tracks);
}

void write_dye_tracks_helper(
        const string& filename,
        unsigned int num_timesteps,
        unsigned int num_channels,
        const vector<SourcedData<DyeTrack, SourceCountHitsList<int>>>&
                dye_tracks) {
    ofstream f(filename);
    f << num_timesteps << "\n";
    f << num_channels << "\n";
    f << dye_tracks.size() << "\n";
    for (const auto& dye_track : dye_tracks) {
        for (unsigned int i = 0; i < num_timesteps * num_channels; i++) {
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

}  // namespace whatprot
