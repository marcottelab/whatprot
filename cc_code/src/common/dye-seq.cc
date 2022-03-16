/******************************************************************************\
* Author: Matthew Beauregard Smith                                             *
* Affiliation: The University of Texas at Austin                               *
* Department: Oden Institute and Institute for Cellular and Molecular Biology  *
* PI: Edward Marcotte                                                          *
* Project: Protein Fluorosequencing                                            *
\******************************************************************************/

// Defining symbols from header:
#include "dye-seq.h"

// Standard C++ library headers:
#include <algorithm>  // needed for std::copy
#include <string>

namespace whatprot {

namespace {
using std::copy;
using std::string;
}  // namespace

DyeSeq::DyeSeq(unsigned int num_channels, const string& s)
        : num_channels(num_channels) {
    length = s.length();
    while (length > 0 && s[length - 1] == '.') {
        length--;
    }
    seq = new short[length];
    for (unsigned int i = 0; i < length; i++) {
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
        : length(other.length), num_channels(other.num_channels) {
    seq = new short[length];
    copy(other.seq, &other.seq[length], seq);
}

DyeSeq::DyeSeq() : seq(NULL), length(0), num_channels(0) {}

DyeSeq::~DyeSeq() {
    if (seq != NULL) {
        delete[] seq;
    }
}

short DyeSeq::operator[](unsigned int i) const {
    if (i < length) {
        return seq[i];
    } else {
        return -1;
    }
}

short& DyeSeq::operator[](unsigned int i) {
    return seq[i];
}

}  // namespace whatprot
