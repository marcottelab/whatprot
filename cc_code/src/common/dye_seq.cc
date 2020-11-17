/******************************************************************************\
* Author: Matthew Beauregard Smith                                             *
* Affiliation: The University of Texas at Austin                               *
* Department: Oden Institute and Institute for Cellular and Molecular Biology  *
* PI: Edward Marcotte                                                          *
* Project: Protein Fluorosequencing                                            *
\******************************************************************************/

// Defining symbols from header:
#include "dye_seq.h"

// Standard C++ library headers:
#include <algorithm>  // needed for std::copy
#include <string>

namespace fluoroseq {

namespace {
using std::copy;
using std::string;
}  // namespace

DyeSeq::DyeSeq(int num_channels, const string& s) : num_channels(num_channels) {
    length = s.length();
    while (s[length - 1] == '.') {
        length--;
    }
    seq = new short[length];
    for (int i = 0; i < length; i++) {
        const char& c = s[i];
        if (c == '.') {
            seq[i] = -1;
        } else if (c >= '0' and c <= '9') {
            // This is a clever trick to convert ASCII code to integer value.
            short sc = (short)(c - '0');
            seq[i] = sc;
        }
    }
}

DyeSeq::DyeSeq(const DyeSeq& other)
        : num_channels(other.num_channels), length(other.length) {
    seq = new short[length];
    copy(other.seq, &other.seq[length], seq);
}

DyeSeq::~DyeSeq() {
    delete[] seq;
}

short DyeSeq::operator[](int i) const {
    if (i < length) {
        return seq[i];
    } else {
        return -1;
    }
}

short& DyeSeq::operator[](int i) {
    return seq[i];
}

}  // namespace fluoroseq
