// Author: Matthew Beauregard Smith (UT Austin)
#include "scored_classification.h"

#include "common/dye_seq.h"

namespace fluoroseq {

ScoredClassification::ScoredClassification(DyeSeq* y, double score)
        : y(y), score(score) {}

}  // namespace fluoroseq
