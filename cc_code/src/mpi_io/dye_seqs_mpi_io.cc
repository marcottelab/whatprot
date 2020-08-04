// Author: Matthew Beauregard Smith (UT Austin)
#include "dye_seqs_mpi_io.h"

#include <fstream>
#include <string>

#include <mpi.h>

#include "common/dye_seq.h"
#include "common/sourced_data.h"

namespace fluoroseq {

namespace {
using std::ifstream;
using std::string;
}  // namespace

void mpi_read_dye_seqs(const std::string& filename,
                       int* num_channels,
                       int* num_dye_seqs,
                       SourcedData<DyeSeq*, SourceCount<int>*>*** dye_seqs) {
    int mpi_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
    if (mpi_rank == 0) {
        mpi_read_dye_seqs_master(filename,
                                 num_channels,
                                 num_dye_seqs,
                                 dye_seqs);
    } else {
        mpi_read_dye_seqs_slave(filename,
                                num_channels,
                                num_dye_seqs,
                                dye_seqs);
    }
}

void mpi_read_dye_seqs_master(
        const std::string& filename,
        int* num_channels,
        int* num_dye_seqs,
        SourcedData<DyeSeq*, SourceCount<int>*>*** dye_seqs) {
    ifstream f(filename);
    f >> *num_channels;
    MPI_Bcast(num_channels,
              1,  // count
              MPI_INT,
              0,  // root
              MPI_COMM_WORLD);
    f >> *num_dye_seqs;
    MPI_Bcast(num_dye_seqs,
              1,  // count
              MPI_INT,
              0,  // root
              MPI_COMM_WORLD);
    string* dye_strings = new string[*num_dye_seqs];
    int* dye_string_lengths = new int[*num_dye_seqs];
    int* dye_seqs_num_peptides = new int[*num_dye_seqs];
    int* dye_seqs_ids = new int[*num_dye_seqs];
    for (int i = 0; i < *num_dye_seqs; i++) {
        f >> dye_strings[i];
        dye_string_lengths[i] = dye_strings[i].length();
        f >> dye_seqs_num_peptides[i];
        f >> dye_seqs_ids[i];
    }
    f.close();
    MPI_Bcast(dye_string_lengths,
              *num_dye_seqs,
              MPI_INT,
              0,  // root
              MPI_COMM_WORLD);
    for (int i = 0; i < *num_dye_seqs; i++) {
        MPI_Bcast(const_cast<char *>(dye_strings[i].c_str()),
                  dye_string_lengths[i],
                  MPI_CHAR,
                  0,  // root
                  MPI_COMM_WORLD);
    }
    delete[] dye_string_lengths;
    MPI_Bcast(dye_seqs_num_peptides,
              *num_dye_seqs,
              MPI_INT,
              0,  // root
              MPI_COMM_WORLD);
    MPI_Bcast(dye_seqs_ids,
              *num_dye_seqs,
              MPI_INT,
              0,  // root
              MPI_COMM_WORLD);
    *dye_seqs = new SourcedData<DyeSeq*, SourceCount<int>*>*[*num_dye_seqs];
    for (int i = 0; i < *num_dye_seqs; i++) {
        (*dye_seqs)[i] = new SourcedData<DyeSeq*, SourceCount<int>*>(
                new DyeSeq(*num_channels, dye_strings[i]),
                new SourceCount<int>(dye_seqs_ids[i],
                                     dye_seqs_num_peptides[i]));
    }
    delete[] dye_strings;
    delete[] dye_seqs_num_peptides;
    delete[] dye_seqs_ids;
}

void mpi_read_dye_seqs_slave(
        const std::string& filename,
        int* num_channels,
        int* num_dye_seqs,
        SourcedData<DyeSeq*, SourceCount<int>*>*** dye_seqs) {
    MPI_Bcast(num_channels,
              1,  // count
              MPI_INT,
              0,  // root
              MPI_COMM_WORLD);
    MPI_Bcast(num_dye_seqs,
              1,  // count
              MPI_INT,
              0,  // root
              MPI_COMM_WORLD);
    int* dye_string_lengths = new int[*num_dye_seqs];
    MPI_Bcast(dye_string_lengths,
              *num_dye_seqs,
              MPI_INT,
              0,  // root
              MPI_COMM_WORLD);
    string* dye_strings = new string[*num_dye_seqs];
    for (int i = 0; i < *num_dye_seqs; i++) {
        char* cstr = new char[dye_string_lengths[i]];
        MPI_Bcast(cstr,
                  dye_string_lengths[i],
                  MPI_CHAR,
                  0,  // root
                  MPI_COMM_WORLD);
        dye_strings[i] = string(cstr, dye_string_lengths[i]);
        delete[] cstr;
    }
    delete[] dye_string_lengths;
    int* dye_seqs_num_peptides = new int[*num_dye_seqs];
    MPI_Bcast(dye_seqs_num_peptides,
              *num_dye_seqs,
              MPI_INT,
              0,  // root
              MPI_COMM_WORLD);
    int* dye_seqs_ids = new int[*num_dye_seqs];
    MPI_Bcast(dye_seqs_ids,
              *num_dye_seqs,
              MPI_INT,
              0,  // root
              MPI_COMM_WORLD);
    *dye_seqs = new SourcedData<DyeSeq*, SourceCount<int>*>*[*num_dye_seqs];
    for (int i = 0; i < *num_dye_seqs; i++) {
        (*dye_seqs)[i] = new SourcedData<DyeSeq*, SourceCount<int>*>(
                new DyeSeq(*num_channels, dye_strings[i]),
                new SourceCount<int>(dye_seqs_ids[i],
                                     dye_seqs_num_peptides[i]));
    }
    delete[] dye_strings;
    delete[] dye_seqs_num_peptides;
    delete[] dye_seqs_ids;
}

}  // namespace fluoroseq
