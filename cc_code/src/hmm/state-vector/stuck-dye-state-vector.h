/******************************************************************************\
* Author: Matthew Beauregard Smith                                             *
* Affiliation: The University of Texas at Austin                               *
* Department: Oden Institute and Institute for Cellular and Molecular Biology  *
* PI: Edward Marcotte                                                          *
* Project: Protein Fluorosequencing                                            *
\******************************************************************************/

#ifndef WHATPROT_HMM_STATE_VECTOR_STUCK_DYE_STATE_VECTOR_H
#define WHATPROT_HMM_STATE_VECTOR_STUCK_DYE_STATE_VECTOR_H

// Standard C++ library headers:
#include <vector>

// Local project headers:
#include "tensor/tensor.h"

namespace whatprot {

class StuckDyeStateVector {
public:
    // To construct a PeptideStateVector, you need to give the order and shape
    // of the underlying tensor.
    StuckDyeStateVector();
    // Put 1.0 in starting state.
    void initialize_from_start();
    // Put 1.0 in every state.
    void initialize_from_finish();
    // Sum of the states.
    double sum() const;
    // Probability of the original source state.
    double source() const;

    double dye;
    double no_dye;
};

}  // namespace whatprot

#endif  // WHATPROT_HMM_STATE_VECTOR_STUCK_DYE_STATE_VECTOR_H
