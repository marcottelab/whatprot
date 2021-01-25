/******************************************************************************\
* Author: Matthew Beauregard Smith                                             *
* Affiliation: The University of Texas at Austin                               *
* Department: Oden Institute and Institute for Cellular and Molecular Biology  *
* PI: Edward Marcotte                                                          *
* Project: Protein Fluorosequencing                                            *
\******************************************************************************/

#ifndef FLUOROSEQ_HMM_STEP_H
#define FLUOROSEQ_HMM_STEP_H

// Local project headers:
#include "tensor/tensor.h"

namespace fluoroseq {

class Step {
public:
    virtual void forward(const Tensor& input,
                         int* edmans,
                         Tensor* output) const = 0;
    virtual void backward(const Tensor& input,
                          int* edmans,
                          Tensor* output) const = 0;
};

}  // namespace fluoroseq

#endif  // FLUOROSEQ_HMM_STEP_H
