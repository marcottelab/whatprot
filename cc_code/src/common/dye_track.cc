// Author: Matthew Beauregard Smith
#include "dye_track.h"
#include "util/vector_hash.h"

#include <algorithm>  // for std::copy
#include <utility>  // for std::move
#include <vector>

namespace fluoroseq {

namespace {
using std::copy;
using std::move;
using std::vector;
}

DyeTrack::DyeTrack(int num_timesteps, int num_channels, const DyeSeq& dye_seq)
        : num_timesteps(num_timesteps), num_channels(num_channels) {
    counts.resize(num_timesteps * num_channels);
    vector<short> cs;
    cs.resize(num_channels);
    for (int t = dye_seq.length - 1; t >= 0; t--) {
        short dye = dye_seq[t];
        if (dye != -1) {
            cs[dye]++;
        }
        if (t < num_timesteps) {
            copy(cs.begin(), cs.end(), &counts[t * num_channels]);
        }
    }
}

DyeTrack::DyeTrack(int num_timesteps, int num_channels, short* counts)
        : num_timesteps(num_timesteps), num_channels(num_channels) {
    this->counts.resize(num_timesteps * num_channels);
    copy(counts, &counts[num_timesteps * num_channels], this->counts.begin());
}

DyeTrack::DyeTrack(int num_timesteps, int num_channels)
        : num_timesteps(num_timesteps), num_channels(num_channels) {
    counts.resize(num_timesteps * num_channels);
}

DyeTrack::DyeTrack(const DyeTrack& other) : num_timesteps(other.num_timesteps),
                                            num_channels(other.num_channels),
                                            counts(other.counts) {}

DyeTrack::DyeTrack(DyeTrack&& other) : num_timesteps(other.num_timesteps),
                                       num_channels(other.num_channels),
                                       counts(move(other.counts)) {}

bool DyeTrack::operator==(const DyeTrack& other) const {
    return num_timesteps == other.num_timesteps
           && num_channels == other.num_channels
           && counts == other.counts;
}

short& DyeTrack::operator()(int t, int c) {
    return counts[t * num_channels + c];
}

short DyeTrack::operator()(int t, int c) const {
    return counts[t * num_channels + c];
}

}  // namespace fluoroseq

namespace std {

size_t hash<fluoroseq::DyeTrack>::operator()(
        const fluoroseq::DyeTrack& dye_track) const {
    return vector_hash(dye_track.counts);
}

}  // namespace std
