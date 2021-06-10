/******************************************************************************\
* Author: Matthew Beauregard Smith                                             *
* Affiliation: The University of Texas at Austin                               *
* Department: Oden Institute and Institute for Cellular and Molecular Biology  *
* PI: Edward Marcotte                                                          *
* Project: Protein Fluorosequencing                                            *
\******************************************************************************/

// Defining symbols from header:
#include "scored-classification.h"

// Standard C++ library headers:
#include <climits>

namespace whatprot {

ScoredClassification::ScoredClassification(int id, double score, double total)
        : score(score), total(total), id(id) {}

ScoredClassification::ScoredClassification()
        : score(INT_MIN), total(0.0), id(-1) {}

double ScoredClassification::adjusted_score() const {
    return score / total;
}

bool operator>(const ScoredClassification& x, const ScoredClassification& y) {
    return (x.score > y.score);
}

}  // namespace whatprot
