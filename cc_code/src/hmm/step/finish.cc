/******************************************************************************\
* Author: Matthew Beauregard Smith                                             *
* Affiliation: The University of Texas at Austin                               *
* Department: Oden Institute and Institute for Cellular and Molecular Biology  *
* PI: Edward Marcotte                                                          *
* Project: Protein Fluorosequencing                                            *
\******************************************************************************/

// Defining symbols from header:
#include "finish.h"

// Local project headers:
#include "hmm/fit/error-model-fitter.h"
#include "tensor/tensor.h"

namespace whatprot {

void Finish::forward(int* edmans, Tensor* output) const {
}

void Finish::backward(const Tensor& input, int* edmans, Tensor* output) const {
    for (int i = 0; i < output->size; i++) {
        output->values[i] = 1.0;
    }
}

void Finish::improve_fit(const Tensor& forward_tensor,
                         const Tensor& backward_tensor,
                         const Tensor& next_backward_tensor,
                         int edmans,
                         double probability,
                         ErrorModelFitter* fitter) const {}

}  // namespace whatprot
