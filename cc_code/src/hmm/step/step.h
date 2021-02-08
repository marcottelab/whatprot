/******************************************************************************\
* Author: Matthew Beauregard Smith                                             *
* Affiliation: The University of Texas at Austin                               *
* Department: Oden Institute and Institute for Cellular and Molecular Biology  *
* PI: Edward Marcotte                                                          *
* Project: Protein Fluorosequencing                                            *
\******************************************************************************/

#ifndef WHATPROT_HMM_STEP_STEP_H
#define WHATPROT_HMM_STEP_STEP_H

// Local project headers:
#include "hmm/fit/error-model-fitter.h"
#include "tensor/tensor.h"

namespace whatprot {

class Step {
public:
    virtual void forward(const Tensor& input,
                         int* edmans,
                         Tensor* output) const = 0;
    virtual void backward(const Tensor& input,
                          int* edmans,
                          Tensor* output) const = 0;
    virtual void improve_fit(const Tensor& forward_tensor,
                             const Tensor& backward_tensor,
                             const Tensor& next_backward_tensor,
                             int edmans,
                             double probability,
                             ErrorModelFitter* fitter) const = 0;
};

}  // namespace whatprot

#endif  // WHATPROT_HMM_STEP_STEP_H
