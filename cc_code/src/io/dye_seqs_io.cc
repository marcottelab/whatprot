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
                   DyeSeq*** dye_seqs,
                   int** dye_seqs_num_peptides,
                   int** dye_seqs_ids) {
    ifstream f(filename);
    f >> *num_channels;
    f >> *num_dye_seqs;
    *dye_seqs = new DyeSeq*[*num_dye_seqs];
    *dye_seqs_num_peptides = new int[*num_dye_seqs];
    *dye_seqs_ids = new int[*num_dye_seqs];
    for (int i = 0; i < *num_dye_seqs; i++) {
        string dye_string;
        f >> dye_string;
        (*dye_seqs)[i] = new DyeSeq(*num_channels, dye_string);
        f >> (*dye_seqs_num_peptides)[i];
        f >> (*dye_seqs_ids)[i];
    }
    f.close();

}  // namespace fluoroseq
