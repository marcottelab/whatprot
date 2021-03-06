/******************************************************************************\
* Author: Matthew Beauregard Smith                                             *
* Affiliation: The University of Texas at Austin                               *
* Department: Oden Institute and Institute for Cellular and Molecular Biology  *
* PI: Edward Marcotte                                                          *
* Project: Protein Fluorosequencing                                            *
\******************************************************************************/

#ifndef WHATPROT_HMM_STATE_VECTOR_PEPTIDE_STATE_VECTOR_H
#define WHATPROT_HMM_STATE_VECTOR_PEPTIDE_STATE_VECTOR_H

// Standard C++ library headers:
#include <vector>

// Local project headers:
#include "tensor/tensor.h"

namespace whatprot {

class PeptideStateVector {
public:
    // To construct a PeptideStateVector, you need to give the order and shape
    // of the underlying tensor.
    PeptideStateVector(int order, const int* shape);
    // Put 1.0 in starting state.
    void initialize_from_start();
    // Put 1.0 in every state.
    void initialize_from_finish();
    // Sum of the states.
    double sum() const;
    // Probability of the original source state.
    double source() const;

    Tensor tensor;
};

}  // namespace whatprot

#endif  // WHATPROT_HMM_STATE_VECTOR_PEPTIDE_STATE_VECTOR_H