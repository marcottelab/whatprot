// Author: Matthew Beauregard Smith (UT Austin)
#ifndef FLUOROSEQ_COMMON_DYE_TRACK_H
#define FLUOROSEQ_COMMON_DYE_TRACK_H

#include "common/dye_seq.h"

namespace fluoroseq {

class DyeTrack {
public:
    DyeTrack(int num_timesteps, int num_channels, const DyeSeq& dye_seq);
    DyeTrack(const DyeTrack& other);
    ~DyeTrack();
    short& operator()(int t, int c);
    short operator()(int t, int c) const;

    short* counts;
    int num_timesteps;
    int num_channels;
};

}  // namespace fluoroseq

#endif  // FLUOROSEQ_COMMON_DYE_TRACK_H
