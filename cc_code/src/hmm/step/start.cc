/******************************************************************************\
* Author: Matthew Beauregard Smith                                             *
* Affiliation: The University of Texas at Austin                               *
* Department: Oden Institute and Institute for Cellular and Molecular Biology  *
* PI: Edward Marcotte                                                          *
* Project: Protein Fluorosequencing                                            *
\******************************************************************************/

// Defining symbols from header:
#include "start.h"

// Local project headers:
#include "hmm/fit/error-model-fitter.h"
#include "tensor/tensor.h"

namespace whatprot {

void Start::forward(const Tensor& input, int* edmans, Tensor* output) const {
    for (int i = 0; i < output->strides[0] - 1; i++) {
        output->values[i] = 0.0;
    }
    output->values[output->strides[0] - 1] = 1.0;
}

void Start::backward(const Tensor& input, int* edmans, Tensor* output) const {}

void Start::improve_fit(const Tensor& forward_tensor,
                        const Tensor& backward_tensor,
                        const Tensor& next_backward_tensor,
                        int edmans,
                        double probability,
                        ErrorModelFitter* fitter) const {}

}  // namespace whatprot
