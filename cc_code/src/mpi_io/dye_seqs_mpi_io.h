// Author: Matthew Beauregard Smith
#ifndef FLUOROSEQ_MPI_IO_DYE_SEQS_MPI_IO_H
#define FLUOROSEQ_MPI_IO_DYE_SEQS_MPI_IO_H

#include <string>

#include "common/dye_seq.h"
#include "common/sourced_data.h"

namespace fluoroseq {

void mpi_read_dye_seqs(const std::string& filename,
                       int* num_channels,
                       int* num_dye_seqs,
                       SourcedData<DyeSeq*, SourceCount<int>*>*** dye_seqs);

void mpi_read_dye_seqs_master(
        const std::string& filename,
        int* num_channels,
        int* num_dye_seqs,
        SourcedData<DyeSeq*, SourceCount<int>*>*** dye_seqs);

void mpi_read_dye_seqs_slave(
        const std::string& filename,
        int* num_channels,
        int* num_dye_seqs,
        SourcedData<DyeSeq*, SourceCount<int>*>*** dye_seqs);

}  // namespace fluoroseq

#endif  // FLUOROSEQ_MPI_IO_DYE_SEQS_MPI_IO_H
