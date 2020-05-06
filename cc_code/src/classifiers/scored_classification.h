// Author: Matthew Beauregard Smith (UT Austin)
#ifndef FLUOROSEQ_CLASSIFIERS_SCORED_CLASSIFICATION_H
#define FLUOROSEQ_CLASSIFIERS_SCORED_CLASSIFICATION_H

#include "common/dye_seq.h"

namespace fluoroseq {

class ScoredClassification {
public:
    ScoredClassification(DyeSeq* y, double score, double total);
    ScoredClassification();
    double adjusted_score();

    DyeSeq* y;  // not owned
    double score;
    double total;
};

}  // namespace fluoroseq

#endif  // FLUOROSEQ_CLASSIFIERS_SCORED_CLASSIFICATION_H
