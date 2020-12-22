/******************************************************************************\
* Author: Matthew Beauregard Smith                                             *
* Affiliation: The University of Texas at Austin                               *
* Department: Oden Institute and Institute for Cellular and Molecular Biology  *
* PI: Edward Marcotte                                                          *
* Project: Protein Fluorosequencing                                            *
\******************************************************************************/

#ifndef FLUOROSEQ_IO_SCORED_CLASSIFICATIONS_IO_H
#define FLUOROSEQ_IO_SCORED_CLASSIFICATIONS_IO_H

// Standard C++ library headers:
#include <string>
#include <vector>

// Local project headers:
#include "common/scored-classification.h"

namespace fluoroseq {

void write_scored_classifications(
        const std::string& filename,
        int total_num_scored_classifications,
        const std::vector<ScoredClassification>& scored_classifications);

void convert_raw_from_scored_classifications(
        const std::vector<ScoredClassification>& scored_classifications,
        int** ids,
        double** scores);

void write_scored_classifications_raw(const std::string& filename,
                                      int num_scored_classifications,
                                      int* ids,
                                      double* scores);

}  // namespace fluoroseq

#endif  // FLUOROSEQ_IO_SCORED_CLASSIFICATIONS_IO_H