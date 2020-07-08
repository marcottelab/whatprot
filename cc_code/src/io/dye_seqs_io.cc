// Author: Matthew Beauregard Smith (UT Austin)
#include "dye_seqs_io.h"

#include <fstream>
#include <string>

#include "common/dye_seq.h"

namespace fluoroseq {

namespace {
using std::ifstream;
using std::string;
}  // namespace

void read_dye_seqs(const string& filename,
                   int* num_channels,
                   int* num_dye_seqs,
                   SourcedData<DyeSeq*, SourceCount<int>*>*** dye_seqs) {
    ifstream f(filename);
    f >> *num_channels;
    f >> *num_dye_seqs;
    *dye_seqs = new SourcedData<DyeSeq*, SourceCount<int>*>*[*num_dye_seqs];
    for (int i = 0; i < *num_dye_seqs; i++) {
        string dye_string;
        f >> dye_string;
        int count;
        f >> count;
        int id;
        f >> id;
        (*dye_seqs)[i] = new SourcedData<DyeSeq*, SourceCount<int>*>(
            new DyeSeq(*num_channels, dye_string),
            new SourceCount<int>(id, count));
    }
    f.close();
}

}  // namespace fluoroseq
