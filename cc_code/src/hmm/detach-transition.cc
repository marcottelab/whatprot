/******************************************************************************\
* Author: Matthew Beauregard Smith                                             *
* Affiliation: The University of Texas at Austin                               *
* Department: Oden Institute and Institute for Cellular and Molecular Biology  *
* PI: Edward Marcotte                                                          *
* Project: Protein Fluorosequencing                                            *
\******************************************************************************/

// Defining symbols from header:
#include "detach-transition.h"

// Local project headers:
#include "tensor/tensor.h"

namespace fluoroseq {

DetachTransition::DetachTransition(double p_detach) : p_detach(p_detach) {}

void DetachTransition::forward(const Tensor& input, int edmans, Tensor* output) const {
    int i_max = (edmans + 1) * input.strides[0];
    double sum = 0.0;
    for (int i = 0; i < i_max; i++) {
        double value = input.values[i];
        output->values[i] = value * (1 - p_detach);
        sum += value;
    }
    output->values[edmans * input.strides[0]] += p_detach * sum;
}

}  // namespace fluoroseq
