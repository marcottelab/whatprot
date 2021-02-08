/******************************************************************************\
* Author: Matthew Beauregard Smith                                             *
* Affiliation: The University of Texas at Austin                               *
* Department: Oden Institute and Institute for Cellular and Molecular Biology  *
* PI: Edward Marcotte                                                          *
* Project: Protein Fluorosequencing                                            *
\******************************************************************************/

// Defining symbols from header:
#include "summation.h"

// Local project headers:
#include "tensor/tensor.h"

namespace whatprot {

Summation::Summation() {}

double Summation::operator()(Tensor* tensor, int timestep) const {
    double sum = 0.0;
    for (int i = 0; i < (timestep + 1) * tensor->strides[0]; i++) {
        sum += tensor->values[i];
    }
    return sum;
}

}  // namespace whatprot
