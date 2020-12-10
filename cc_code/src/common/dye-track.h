/******************************************************************************\
* Author: Matthew Beauregard Smith                                             *
* Affiliation: The University of Texas at Austin                               *
* Department: Oden Institute and Institute for Cellular and Molecular Biology  *
* PI: Edward Marcotte                                                          *
* Project: Protein Fluorosequencing                                            *
\******************************************************************************/

#ifndef FLUOROSEQ_COMMON_DYE_TRACK_H
#define FLUOROSEQ_COMMON_DYE_TRACK_H

// Standard C++ library headers:
#include <functional>
#include <vector>

// Local project headers:
#include "common/dye-seq.h"
#include "util/vector-hash.h"

namespace fluoroseq {

class DyeTrack {
public:
    DyeTrack(int num_timesteps, int num_channels, const DyeSeq& dye_seq);
    DyeTrack(int num_timesteps, int num_channels, short* counts);
    DyeTrack(int num_timesteps, int num_channels);
    DyeTrack(const DyeTrack& other);
    DyeTrack(DyeTrack&& other);
    DyeTrack& operator=(DyeTrack&& other);
    bool operator==(const DyeTrack& other) const;
    short& operator()(int t, int c);
    short operator()(int t, int c) const;

    std::vector<short> counts;
    int num_timesteps;
    int num_channels;
};

}  // namespace fluoroseq

// An std namespace injection is the accepted way of creating a new hash
// function.
namespace std {

template <>
struct hash<fluoroseq::DyeTrack> {
public:
    size_t operator()(const fluoroseq::DyeTrack& dye_track) const;

    hash<vector<short>> vector_hash;
};

}  // namespace std

#endif  // FLUOROSEQ_COMMON_DYE_TRACK_H
