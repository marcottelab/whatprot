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

PeptideStateVector::PeptideStateVector(unsigned int order,
                                       const unsigned int* shape)
        : tensor(order, shape), p_detached(0.0), allow_detached(true) {
    for (unsigned int o = 0; o < order; o++) {
        range.min.push_back(0);
        range.max.push_back(shape[o]);
    }
}

PeptideStateVector::PeptideStateVector(const KDRange& range)
        : tensor(range), p_detached(0.0), allow_detached(true) {}

void PeptideStateVector::initialize_from_start() {
    tensor.values[tensor.strides[0] - 1] = 1.0;
}

void PeptideStateVector::initialize_from_finish() {
    for (unsigned int i = 0; i < tensor.size; i++) {
        tensor.values[i] = 1.0;
    }
    p_detached = 1.0;
}

double PeptideStateVector::sum() const {
    return tensor.sum(range) + p_detached;
}

double PeptideStateVector::source() const {
    return tensor.values[tensor.strides[0] - 1];
}

}  // namespace whatprot
