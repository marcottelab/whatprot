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

void Start::forward(int* edmans, Tensor* tsr) const {
    tsr->values[tsr->strides[0] - 1] = 1.0;
}

void Start::backward(const Tensor& input, int* edmans, Tensor* output) const {
    if (&input != output) {
        for (int i = 0; i < output->size; i++) {
            output->values[i] = input.values[i];
        }
    }
}

void Start::improve_fit(const Tensor& forward_tensor,
                        const Tensor& backward_tensor,
                        const Tensor& next_backward_tensor,
                        int edmans,
                        double probability,
                        ErrorModelFitter* fitter) const {}

}  // namespace whatprot
