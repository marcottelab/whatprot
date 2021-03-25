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

void Finish::forward(int* num_edmans, PeptideStateVector* psv) const {}

void Finish::backward(const PeptideStateVector& input, int* num_edmans,
                      PeptideStateVector* output) const {
    for (int i = 0; i < output->tensor.size; i++) {
        output->tensor.values[i] = 1.0;
    }
    // The zeroth dimension in the Tensors is time, so output->shape[0] gives
    // the total number of timesteps being considered. There is one less Edman
    // than the number of timesteps, because no Edman is done before the zeroth
    // timestep.
    (*num_edmans) = output->tensor.shape[0] - 1;
}

void Finish::improve_fit(const PeptideStateVector& forward_psv,
                         const PeptideStateVector& backward_psv,
                         const PeptideStateVector& next_backward_psv, int num_edmans,
                         double probability,
                         ErrorModelFitter* fitter) const {}

}  // namespace whatprot
