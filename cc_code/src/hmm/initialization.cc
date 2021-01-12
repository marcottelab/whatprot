/******************************************************************************\
* Author: Matthew Beauregard Smith                                             *
* Affiliation: The University of Texas at Austin                               *
* Department: Oden Institute and Institute for Cellular and Molecular Biology  *
* PI: Edward Marcotte                                                          *
* Project: Protein Fluorosequencing                                            *
\******************************************************************************/

// Defining symbols from header:
#include "initialization.h"

// Local project headers:
#include "tensor/tensor.h"

namespace fluoroseq {

void Initialization::forward(Tensor* tensor) const {
    for (int i = 0; i < tensor->strides[0] - 1; i++) {
        tensor->values[i] = 0.0;
    }
    tensor->values[tensor->strides[0] - 1] = 1.0;
}

}  // namespace fluoroseq
