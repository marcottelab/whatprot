/******************************************************************************\
* Author: Matthew Beauregard Smith                                             *
* Affiliation: The University of Texas at Austin                               *
* Department: Oden Institute and Institute for Cellular and Molecular Biology  *
* PI: Edward Marcotte                                                          *
* Project: Protein Fluorosequencing                                            *
\******************************************************************************/

// You should only be including this for MPI enabled builds.
#ifndef FLUOROSEQ_IO_MPI_COUNTS_DISPLS_H
#define FLUOROSEQ_IO_MPI_COUNTS_DISPLS_H

namespace fluoroseq {

#ifdef USE_MPI
void mpi_counts_displs(int total_count,
                       int block_size,
                       int** counts,
                       int** displs);
#endif  // USE_MPI

}  // namespace fluoroseq

#endif  // FLUOROSEQ_IO_MPI_COUNTS_DISPLS_H