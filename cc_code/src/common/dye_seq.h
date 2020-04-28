// Author: Matthew Beauregard Smith (UT Austin)
#ifndef FLUOROSEQ_COMMON_DYE_SEQ_H
#define FLUOROSEQ_COMMON_DYE_SEQ_H

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
    DyeSeq(int num_channels, const string& s, int num_peptides, int id);
    DyeSeq(const DyeSeq& other);
    ~DyeSeq();
    short operator[](int i) const;

    short* seq;
    int length;
    int num_channels;
    int num_peptides;  // number of peptides that map to this dye seq.
    int id;  // a unique ID for this dye seq.
};

}  // namespace fluoroseq

#endif  // FLUOROSEQ_COMMON_DYE_SEQ_H
