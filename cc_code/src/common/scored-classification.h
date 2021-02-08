/******************************************************************************\
* Author: Matthew Beauregard Smith                                             *
* Affiliation: The University of Texas at Austin                               *
* Department: Oden Institute and Institute for Cellular and Molecular Biology  *
* PI: Edward Marcotte                                                          *
* Project: Protein Fluorosequencing                                            *
\******************************************************************************/

#ifndef WHATPROT_COMMON_SCORED_CLASSIFICATION_H
#define WHATPROT_COMMON_SCORED_CLASSIFICATION_H

namespace whatprot {

class ScoredClassification {
public:
    ScoredClassification(int id, double score, double total);
    ScoredClassification();
    double adjusted_score() const;

    double score;
    double total;
    int id;
};

bool operator>(const ScoredClassification& x, const ScoredClassification& y);

}  // namespace whatprot

#endif  // WHATPROT_COMMON_SCORED_CLASSIFICATION_H
