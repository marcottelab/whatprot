// Author: Matthew Beauregard Smith (UT Austin)
//
// For MPI version, define compiler macro USE_MPI when building.
#include "dye_seqs_io.h"

#include <algorithm>
#include <fstream>
#include <string>
#include <vector>

#ifdef USE_MPI
#include <mpi.h>
#endif  // USE_MPI

#include "common/dye_seq.h"
#include "common/sourced_data.h"

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
                   vector<SourcedData<DyeSeq, SourceCount<int>>>* dye_seqs,
                   MpiReadMode mode) {
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
#ifdef USE_MPI
    switch (mode) {
        case BROADCAST:
            broadcast_dye_seqs(num_channels,
                               total_num_dye_seqs,
                               &dye_string_lengths,
                               &dye_strings,
                               &dye_seqs_num_peptides,
                               &dye_seqs_ids);
            num_dye_seqs = *total_num_dye_seqs;
        break;
        case SCATTER:
            scatter_dye_seqs(num_channels,
                             total_num_dye_seqs,
                             &num_dye_seqs,
                             &dye_string_lengths,
                             &dye_strings,
                             &dye_seqs_num_peptides,
                             &dye_seqs_ids);
        break;
    }
#endif  // USE_MPI
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
#ifdef USE_MPI
    int mpi_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
    if (mpi_rank != 0) {
        return;
    }
#endif  // USE_MPI
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

#ifdef USE_MPI

void broadcast_dye_seqs(int* num_channels,
                        int* num_dye_seqs,
                        int** dye_string_lengths,
                        char*** dye_strings,
                        int** dye_seqs_num_peptides,
                        int** dye_seqs_ids) {
    int mpi_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
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
    if (mpi_rank != 0) {
        *dye_string_lengths = new int[*num_dye_seqs];
        *dye_strings = new char*[*num_dye_seqs];
        *dye_seqs_num_peptides = new int[*num_dye_seqs];
        *dye_seqs_ids = new int[*num_dye_seqs];
    }
    MPI_Bcast(*dye_string_lengths,
              *num_dye_seqs,
              MPI_INT,
              0,  // root
              MPI_COMM_WORLD);
    for (int i = 0; i < *num_dye_seqs; i++) {
        if (mpi_rank != 0) {
            (*dye_strings)[i] = new char[(*dye_string_lengths)[i]];
        }
        MPI_Bcast((*dye_strings)[i],
                  (*dye_string_lengths)[i],
                  MPI_CHAR,
                  0,  // root
                  MPI_COMM_WORLD);
    }
    MPI_Bcast(*dye_seqs_num_peptides,
              *num_dye_seqs,
              MPI_INT,
              0,  // root
              MPI_COMM_WORLD);
    MPI_Bcast(*dye_seqs_ids,
              *num_dye_seqs,
              MPI_INT,
              0,  // root
              MPI_COMM_WORLD);
}

void scatter_dye_seqs(int* num_channels,
                      int* total_num_dye_seqs,
                      int* num_dye_seqs,
                      int** dye_string_lengths,
                      char*** dye_strings,
                      int** dye_seqs_num_peptides,
                      int** dye_seqs_ids) {
    int mpi_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
    int mpi_size;
    MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);
    MPI_Bcast(num_channels,
              1,  // count
              MPI_INT,
              0,  // root
              MPI_COMM_WORLD);
    MPI_Bcast(total_num_dye_seqs,
              1,  // count
              MPI_INT,
              0,  // root
              MPI_COMM_WORLD);
    int* mpi_counts;
    int* mpi_displs;
    mpi_counts_displs(*total_num_dye_seqs,
                      1,  // block size
                      &mpi_counts,
                      &mpi_displs);
    *num_dye_seqs = mpi_counts[mpi_rank];
    int* string_lengths_recv = new int[*num_dye_seqs];
    if (mpi_rank != 0) {
        *dye_strings = new char*[*num_dye_seqs];
    }
    // if (mpi_rank != 0) {
    //     *dye_string_lengths = new int[*num_dye_seqs];
    //     *dye_strings = new char*[*num_dye_seqs];
    //     *dye_seqs_num_peptides = new int[*num_dye_seqs];
    //     *dye_seqs_ids = new int[*num_dye_seqs];
    // }
    MPI_Scatterv(*dye_string_lengths,  // sendbuf
                 mpi_counts,  // sendcounts
                 mpi_displs,  // displs
                 MPI_INT,  // sendtype
                 string_lengths_recv,  // recvbuf
                 *num_dye_seqs,  // recvcount
                 MPI_INT,  // recvtype
                 0,  // root
                 MPI_COMM_WORLD);
    if (mpi_rank == 0) {
        delete[] *dye_string_lengths;
    }
    *dye_string_lengths = string_lengths_recv;
    if (mpi_rank == 0) {
        for (int r = 1; r < mpi_size; r++) {
            for (int i = 0; i < *num_dye_seqs; i++) {
                int index = mpi_displs[r] + i;
                MPI_Send((*dye_strings)[index],  // buf
                         (*dye_string_lengths)[index],  // count
                         MPI_CHAR,
                         r,  // dest
                         0,  // tag
                         MPI_COMM_WORLD);
                delete[] (*dye_strings)[index];
            }
        }
    } else {
        for (int i = 0; i < *num_dye_seqs; i++) {
            (*dye_strings)[i] = new char[(*dye_string_lengths)[i]];
            MPI_Recv((*dye_strings)[i],  // buf
                     (*dye_string_lengths)[i],  // count
                     MPI_CHAR,
                     0,  // source
                     0,  // tag
                     MPI_COMM_WORLD);
        }
    }
    int* num_peptides_recv = new int[*num_dye_seqs];
    MPI_Scatterv(*dye_seqs_num_peptides,  // sendbuf
                 mpi_counts,  // sendcounts
                 mpi_displs,  // displs
                 MPI_INT,  // sendtype
                 num_peptides_recv,  // recvbuf
                 *num_dye_seqs,  // recvcount
                 MPI_INT,  // recvtype
                 0,  // root
                 MPI_COMM_WORLD);
    if (mpi_rank == 0) {
        delete[] *dye_seqs_num_peptides;
    }
    *dye_seqs_num_peptides = num_peptides_recv;
    int* ids_recv = new int[*num_dye_seqs];
    MPI_Scatterv(*dye_seqs_ids,  // sendbuf
                 mpi_counts,  // sendcounts
                 mpi_displs,  // displs
                 MPI_INT,  // sendtype
                 ids_recv,  // recvbuf
                 *num_dye_seqs,  // recvcount
                 MPI_INT,  // recvtype
                 0,  // root
                 MPI_COMM_WORLD);
    if (mpi_rank == 0) {
        delete[] *dye_seqs_ids;
    }
    *dye_seqs_ids = ids_recv;
    delete[] mpi_counts;
    delete[] mpi_displs;
}

#endif  // USE_MPI

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
        dye_seqs->push_back(
                SourcedData<DyeSeq, SourceCount<int>>(
                        DyeSeq(num_channels, dye_string),
                        SourceCount<int>(dye_seqs_ids[i],
                                         dye_seqs_num_peptides[i])));
    }
}

}  // namespace fluoroseq