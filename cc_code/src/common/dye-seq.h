/******************************************************************************\
* Author: Matthew Beauregard Smith                                             *
* Affiliation: The University of Texas at Austin                               *
* Department: Oden Institute and Institute for Cellular and Molecular Biology  *
* PI: Edward Marcotte                                                          *
* Project: Protein Fluorosequencing                                            *
\******************************************************************************/

#ifndef WHATPROT_COMMON_DYE_SEQ_H
#define WHATPROT_COMMON_DYE_SEQ_H

// Standard C++ library headers:
#include <string>

namespace whatprot {

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
    DyeSeq(unsigned int num_channels, const string& s);
    DyeSeq(const DyeSeq& other);
    ~DyeSeq();
    short operator[](unsigned int i) const;
    short& operator[](unsigned int i);

    short* seq;
    unsigned int length;
    unsigned int num_channels;
};

}  // namespace whatprot

#endif  // WHATPROT_COMMON_DYE_SEQ_H
