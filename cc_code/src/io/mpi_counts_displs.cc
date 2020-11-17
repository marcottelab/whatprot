/******************************************************************************\
* Author: Matthew Beauregard Smith                                             *
* Affiliation: The University of Texas at Austin                               *
* Department: Oden Institute and Institute for Cellular and Molecular Biology  *
* PI: Edward Marcotte                                                          *
* Project: Protein Fluorosequencing                                            *
\******************************************************************************/

// You should only be including this for MPI enabled builds.

// Defining symbols from header:
#include "mpi_counts_displs.h"

// MPI header:
#ifdef USE_MPI
#include <mpi.h>
#endif  // USE_MPI

namespace fluoroseq {

#ifdef USE_MPI
void mpi_counts_displs(int total_count,
                       int block_size,
                       int** counts,
                       int** displs) {
    int mpi_size;
    MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);
    *counts = new int[mpi_size];
    *displs = new int[mpi_size];
    for (int i = 0; i < mpi_size; i++) {
        int begin = (long)total_count * (long)i / (long)mpi_size;
        int end = (long)total_count * (long)(i + 1) / (long)mpi_size;
        (*counts)[i] = (end - begin) * block_size;
        (*displs)[i] = begin * block_size;
    }
}
#endif  // USE_MPI

}  // namespace fluoroseq
