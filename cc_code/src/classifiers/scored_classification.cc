// Author: Matthew Beauregard Smith (UT Austin)
#include "scored_classification.h"

#include <climits>

#include "common/dye_seq.h"

namespace fluoroseq {

ScoredClassification::ScoredClassification(DyeSeq* y,
                                           double score,
                                           double total) : y(y),
                                                           score(score),
                                                           total(total) {}

ScoredClassification::ScoredClassification() : y(NULL),
                                               score(INT_MIN),
                                               total(0.0) {}

double ScoredClassification::adjusted_score() {
    return score / total;
}

ScoredClassification* merge_scores(ScoredClassification* a,
                                   ScoredClassification* b) {
    if (a->score > b->score) {
        a->total += b->total;
        delete b;
        return a;
    } else {
        b->total += a->total;
        delete a;
        return b;
    }
}

}  // namespace fluoroseq
