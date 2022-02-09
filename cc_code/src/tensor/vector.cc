/******************************************************************************\
* Author: Matthew Beauregard Smith                                             *
* Affiliation: The University of Texas at Austin                               *
* Department: Oden Institute and Institute for Cellular and Molecular Biology  *
* PI: Edward Marcotte                                                          *
* Project: Protein Fluorosequencing                                            *
\******************************************************************************/

// File under test:
#include "vector.h"

namespace whatprot {

Vector::Vector(unsigned int min, unsigned int max, int stride, double* values)
        : values(values), min(min), max(max), stride(stride) {}

double& Vector::operator[](int i) {
    return values[(i - min) * stride];
}

double Vector::operator[](int i) const {
    return values[(i - min) * stride];
}

}  // namespace whatprot
