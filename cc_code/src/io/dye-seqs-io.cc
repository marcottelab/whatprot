/******************************************************************************\
* Author: Matthew Beauregard Smith                                             *
* Affiliation: The University of Texas at Austin                               *
* Department: Oden Institute and Institute for Cellular and Molecular Biology  *
* PI: Edward Marcotte                                                          *
* Project: Protein Fluorosequencing                                            *
\******************************************************************************/

// Defining symbols from header:
#include "dye-seqs-io.h"

// Standard C++ library headers:
#include <algorithm>
#include <fstream>
#include <string>
#include <vector>

// Local project headers:
#include "common/dye-seq.h"
#include "common/sourced-data.h"

namespace fluoroseq {

namespace {
using std::copy;
using std::ifstream;
using std::string;
using std::vector;
}  // namespace

void read_dye_seqs(const string& filename,
                   int* num_channels,
                   int* total_num_dye_seqs,
                   vector<SourcedData<DyeSeq, SourceCount<int>>>* dye_seqs) {
    int* dye_string_lengths;
    char** dye_strings;
    int* dye_seqs_num_peptides;
    int* dye_seqs_ids;
    read_dye_seqs_raw(filename,
                      num_channels,
                      total_num_dye_seqs,
                      &dye_string_lengths,
                      &dye_strings,
                      &dye_seqs_num_peptides,
                      &dye_seqs_ids);
    int num_dye_seqs = *total_num_dye_seqs;
    convert_dye_seqs_from_raw(*num_channels,
                              num_dye_seqs,
                              dye_string_lengths,
                              dye_strings,
                              dye_seqs_num_peptides,
                              dye_seqs_ids,
                              dye_seqs);
    delete[] dye_string_lengths;
    for (int i = 0; i < num_dye_seqs; i++) {
        delete[] dye_strings[i];
    }
    delete[] dye_strings;
    delete[] dye_seqs_num_peptides;
    delete[] dye_seqs_ids;
}

void read_dye_seqs_raw(const string& filename,
                       int* num_channels,
                       int* num_dye_seqs,
                       int** dye_string_lengths,
                       char*** dye_strings,
                       int** dye_seqs_num_peptides,
                       int** dye_seqs_ids) {
    ifstream f(filename);
    f >> *num_channels;
    f >> *num_dye_seqs;
    *dye_string_lengths = new int[*num_dye_seqs];
    *dye_strings = new char*[*num_dye_seqs];
    *dye_seqs_num_peptides = new int[*num_dye_seqs];
    *dye_seqs_ids = new int[*num_dye_seqs];
    for (int i = 0; i < *num_dye_seqs; i++) {
        string dye_string;
        f >> dye_string;
        (*dye_string_lengths)[i] = dye_string.length();
        (*dye_strings)[i] = new char[dye_string.length()];
        copy(dye_string.c_str(),
             dye_string.c_str() + dye_string.length(),
             (*dye_strings)[i]);
        f >> (*dye_seqs_num_peptides)[i];
        f >> (*dye_seqs_ids)[i];
    }
    f.close();
}

void convert_dye_seqs_from_raw(
        int num_channels,
        int num_dye_seqs,
        int* dye_string_lengths,
        char** dye_strings,
        int* dye_seqs_num_peptides,
        int* dye_seqs_ids,
        vector<SourcedData<DyeSeq, SourceCount<int>>>* dye_seqs) {
    dye_seqs->reserve(num_dye_seqs);
    for (int i = 0; i < num_dye_seqs; i++) {
        string dye_string(dye_strings[i], dye_string_lengths[i]);
        dye_seqs->push_back(SourcedData<DyeSeq, SourceCount<int>>(
                DyeSeq(num_channels, dye_string),
                SourceCount<int>(dye_seqs_ids[i], dye_seqs_num_peptides[i])));
    }
}

}  // namespace fluoroseq