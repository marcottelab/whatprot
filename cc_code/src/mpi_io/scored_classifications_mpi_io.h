// Author: Matthew Beauregard Smith (UT Austin)
#ifndef FLUOROSEQ_MPI_IO_SCORED_CLASSIFICATIONS_MPI_IO
#define FLUOROSEQ_MPI_IO_SCORED_CLASSIFICATIONS_MPI_IO

#include <string>

#include "common/scored_classification.h"

namespace fluoroseq {

void mpi_write_scored_classifications(
        const std::string& filename,
        int total_num_scored_classifications,
        int num_scored_classifications,
        const ScoredClassification* scored_classifications);

void mpi_write_scored_classifications_master(
        const std::string& filename,
        int total_num_scored_classifications,
        int num_scored_classifications,
        const ScoredClassification* scored_classifications);

void mpi_write_scored_classifications_slave(
        const std::string& filename,
        int total_num_scored_classifications,
        int num_scored_classifications,
        const ScoredClassification* scored_classifications);

}  // namespace fluoroseq

#endif  // FLUOROSEQ_MPI_IO_SCORED_CLASSIFICATIONS_MPI_IO