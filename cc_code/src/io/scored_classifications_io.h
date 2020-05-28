// Author: Matthew Beauregard Smith (UT Austin)
#ifndef FLUOROSEQ_IO_SCORED_CLASSIFICATIONS_IO_H
#define FLUOROSEQ_IO_SCORED_CLASSIFICATIONS_IO_H

#include <string>

#include "common/scored_classification.h"

namespace fluoroseq {

void write_scored_classifications(
        const std::string& filename,
        int num_scored_classifications,
        const ScoredClassification* scored_classifications);

}  // namespace fluoroseq

#endif  // FLUOROSEQ_IO_SCORED_CLASSIFICATIONS_IO_H