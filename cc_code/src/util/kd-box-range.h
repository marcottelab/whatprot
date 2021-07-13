/******************************************************************************\
* Author: Matthew Beauregard Smith                                             *
* Affiliation: The University of Texas at Austin                               *
* Department: Oden Institute and Institute for Cellular and Molecular Biology  *
* PI: Edward Marcotte                                                          *
* Project: Protein Fluorosequencing                                            *
\******************************************************************************/

#ifndef WHATPROT_UTIL_KD_BOX_RANGE_H
#define WHATPROT_UTIL_KD_BOX_RANGE_H

// Standard C++ library headers:
#include <vector>

namespace whatprot {

class KDBoxRange {
public:
    KDBoxRange* intersect(const KDBoxRange& other) const;
    bool is_empty() const;

    std::vector<unsigned int> min;
    std::vector<unsigned int> max;
};

}  // namespace whatprot

#endif  // WHATPROT_UTIL_KD_BOX_RANGE_H
