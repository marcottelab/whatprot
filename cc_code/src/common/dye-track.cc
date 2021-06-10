/******************************************************************************\
* Author: Matthew Beauregard Smith                                             *
* Affiliation: The University of Texas at Austin                               *
* Department: Oden Institute and Institute for Cellular and Molecular Biology  *
* PI: Edward Marcotte                                                          *
* Project: Protein Fluorosequencing                                            *
\******************************************************************************/

// Defining symbols from header:
#include "dye-track.h"
#include "util/vector-hash.h"

// Standard C++ library headers:
#include <algorithm>  // for std::copy
#include <utility>  // for std::move
#include <vector>

namespace whatprot {

namespace {
using std::copy;
using std::move;
using std::vector;
}  // namespace

DyeTrack::DyeTrack(unsigned int num_timesteps, unsigned int num_channels, const DyeSeq& dye_seq)
        : num_timesteps(num_timesteps), num_channels(num_channels) {
    counts.resize(num_timesteps * num_channels);
    vector<short> cs;
    cs.resize(num_channels);
    // Need t as signed int because we decrement it towards 0.
    for (int t = dye_seq.length - 1; t >= 0; t--) {
        short dye = dye_seq[t];
        if (dye != -1) {
            cs[dye]++;
        }
        if (t < (int) num_timesteps) {
            copy(cs.begin(), cs.end(), &counts[t * num_channels]);
        }
    }
}

DyeTrack::DyeTrack(unsigned int num_timesteps, unsigned int num_channels, short* counts)
        : num_timesteps(num_timesteps), num_channels(num_channels) {
    this->counts.resize(num_timesteps * num_channels);
    copy(counts, &counts[num_timesteps * num_channels], this->counts.begin());
}

DyeTrack::DyeTrack(unsigned int num_timesteps, unsigned int num_channels)
        : num_timesteps(num_timesteps), num_channels(num_channels) {
    counts.resize(num_timesteps * num_channels);
}

DyeTrack::DyeTrack(const DyeTrack& other)
        : counts(other.counts),
          num_timesteps(other.num_timesteps),
          num_channels(other.num_channels) {}

DyeTrack::DyeTrack(DyeTrack&& other)
        : counts(move(other.counts)),
          num_timesteps(other.num_timesteps),
          num_channels(other.num_channels) {}

DyeTrack& DyeTrack::operator=(DyeTrack&& other) {
    num_timesteps = other.num_timesteps;
    num_channels = other.num_channels;
    counts = move(other.counts);
    return *this;
}

bool DyeTrack::operator==(const DyeTrack& other) const {
    return num_timesteps == other.num_timesteps
           && num_channels == other.num_channels && counts == other.counts;
}

short& DyeTrack::operator()(int t, int c) {
    return counts[t * num_channels + c];
}

short DyeTrack::operator()(int t, int c) const {
    return counts[t * num_channels + c];
}

}  // namespace whatprot

namespace std {

size_t hash<whatprot::DyeTrack>::operator()(
        const whatprot::DyeTrack& dye_track) const {
    return vector_hash(dye_track.counts);
}

}  // namespace std
