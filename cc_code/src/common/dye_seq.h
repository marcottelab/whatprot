/******************************************************************************\
* Author: Matthew Beauregard Smith                                             *
* Affiliation: The University of Texas at Austin                               *
* Department: Oden Institute and Institute for Cellular and Molecular Biology  *
* PI: Edward Marcotte                                                          *
* Project: Protein Fluorosequencing                                            *
\******************************************************************************/

#ifndef FLUOROSEQ_COMMON_DYE_SEQ_H
#define FLUOROSEQ_COMMON_DYE_SEQ_H

// Standard C++ library headers:
#include <string>

namespace fluoroseq {

namespace {
using std::string;
}  // namespace

class DyeSeq {
  public:
    // Parameters:
    //     -- num_channels is a 1 digit positive non-zero integer.
    //     -- s is a string representation of the dye sequence, and should
    //        consist only of '.' and character representations of integers that
    //        are strictly less than num_channels.
    DyeSeq(int num_channels, const string& s);
    DyeSeq(const DyeSeq& other);
    ~DyeSeq();
    short operator[](int i) const;
    short& operator[](int i);

    short* seq;
    int length;
    int num_channels;
};

}  // namespace fluoroseq

#endif  // FLUOROSEQ_COMMON_DYE_SEQ_H
