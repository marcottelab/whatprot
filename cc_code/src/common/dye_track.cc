// Author: Matthew Beauregard Smith
#include "dye_track.h"

#include <algorithm>  // for std::copy

namespace fluoroseq {

namespace {
using std::copy;
}

DyeTrack::DyeTrack(int num_timesteps, int num_channels, const DyeSeq& dye_seq)
        : num_timesteps(num_timesteps), num_channels(num_channels) {
    counts = new short[num_timesteps * num_channels]();
    short* cs = new short[num_channels]();
    for (int t = dye_seq.length - 1; t >= 0; t--) {
        short dye = dye_seq[t];
        if (dye != -1) {
            cs[dye]++;
        }
        if (t < num_timesteps) {
            copy(cs, &cs[num_channels], &counts[t * num_channels]);
        }
    }
    delete[] cs;
}

DyeTrack::DyeTrack(int num_timesteps, int num_channels)
        : num_timesteps(num_timesteps), num_channels(num_channels) {
    counts = new short[num_timesteps * num_channels]();
}

DyeTrack::DyeTrack(const DyeTrack& other) : num_timesteps(other.num_timesteps),
                                            num_channels(other.num_channels) {
    counts = new short[num_timesteps * num_channels];
    copy(other.counts, &other.counts[num_timesteps * num_channels], counts);
}

DyeTrack::~DyeTrack() {
    delete[] counts;
}

short& DyeTrack::operator()(int t, int c) {
    return counts[t * num_channels + c];
}

short DyeTrack::operator()(int t, int c) const {
    return counts[t * num_channels + c];
}

}  // namespace fluoroseq
