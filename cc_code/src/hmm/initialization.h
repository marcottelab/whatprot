/******************************************************************************\
* Author: Matthew Beauregard Smith                                             *
* Affiliation: The University of Texas at Austin                               *
* Department: Oden Institute and Institute for Cellular and Molecular Biology  *
* PI: Edward Marcotte                                                          *
* Project: Protein Fluorosequencing                                            *
\******************************************************************************/

#ifndef FLUOROSEQ_HMM_INITIALIZATION_H
#define FLUOROSEQ_HMM_INITIALIZATION_H

// Local project headers:
#include "tensor/tensor.h"

namespace fluoroseq {

class Initialization {
public:
    void forward(Tensor* tensor) const;
};

}  // namespace fluoroseq

#endif  // FLUOROSEQ_HMM_INITIALIZATION_H