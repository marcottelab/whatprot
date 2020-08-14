// Author: Matthew Beauregard Smith (UT Austin)
#ifndef FLUOROSEQ_IO_SCORED_CLASSIFICATIONS_IO_H
#define FLUOROSEQ_IO_SCORED_CLASSIFICATIONS_IO_H

#include <string>

#include "common/scored_classification.h"

namespace fluoroseq {

void write_scored_classifications(
        const std::string& filename,
        int total_num_scored_classifications,
        int num_scored_classifications,
        const ScoredClassification* scored_classifications);

void convert_raw_from_scored_classifications(
        int num_scored_classifications,
        const ScoredClassification* scored_classifications,
        int** ids,
        double** scores);

#ifdef USE_MPI
void gather_scored_classifications(int total_num_scored_classifications,
                                   int num_scored_classifications,
                                   int** ids,
                                   double** scores);
#endif  // USE_MPI

void write_scored_classifications_raw(const std::string& filename,
                                      int num_scored_classifications,
                                      int* ids,
                                      double* scores);

}  // namespace fluoroseq

#endif  // FLUOROSEQ_IO_SCORED_CLASSIFICATIONS_IO_H