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

void Start::forward(int* num_edmans, PeptideStateVector* psv) const {
    psv->tensor.values[psv->tensor.strides[0] - 1] = 1.0;
}

void Start::backward(const PeptideStateVector& input, int* num_edmans,
                     PeptideStateVector* output) const {
    if (&input != output) {
        for (int i = 0; i < output->tensor.size; i++) {
            output->tensor.values[i] = input.tensor.values[i];
        }
    }
    
}

void Start::improve_fit(const PeptideStateVector& forward_psv,
                        const PeptideStateVector& backward_psv,
                        const PeptideStateVector& next_backward_psv, int num_edmans,
                        double probability,
                        ErrorModelFitter* fitter) const {}

}  // namespace whatprot
