/******************************************************************************\
* Author: Matthew Beauregard Smith                                             *
* Affiliation: The University of Texas at Austin                               *
* Department: Oden Institute and Institute for Cellular and Molecular Biology  *
* PI: Edward Marcotte                                                          *
* Project: Protein Fluorosequencing                                            *
\******************************************************************************/

// Defining symbols from header:
#include "peptide-state-vector.h"

// Standard C++ library headers:
#include <vector>

// Local project headers:
#include "tensor/tensor.h"

namespace whatprot {

PeptideStateVector::PeptideStateVector(int order, const int* shape)
        : tensor(order, shape) {}

void PeptideStateVector::initialize_from_start() {
    tensor.values[tensor.strides[0] - 1] = 1.0;
}

void PeptideStateVector::initialize_from_finish() {
    for (int i = 0; i < tensor.size; i++) {
        tensor.values[i] = 1.0;
    }
}

double PeptideStateVector::sum() const {
    return tensor.sum();
}

double PeptideStateVector::source() const {
    return tensor.values[tensor.strides[0] - 1];
}

}  // namespace whatprot
