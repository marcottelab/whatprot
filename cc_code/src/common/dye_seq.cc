// Author: Matthew Beauregard Smith (UT Austin)
#include "dye_seq.h"

#include <algorithm>  // needed for std::copy
#include <string>

namespace fluoroseq {

namespace {
using std::copy;
using std::string;
}  // namespace

DyeSeq::DyeSeq(int num_channels, const string& s, int num_peptides, int id)
        : num_channels(num_channels), num_peptides(num_peptides), id(id) {
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
            short sc = (short) (c - '0');
            seq[i] = sc;
        }
    }
}

DyeSeq::DyeSeq(const DyeSeq& other) : num_channels(other.num_channels),
                                      length(other.length) {
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

}  // namespace fluoroseq
