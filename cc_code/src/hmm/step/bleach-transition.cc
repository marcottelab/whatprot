/******************************************************************************\
* Author: Matthew Beauregard Smith                                             *
* Affiliation: The University of Texas at Austin                               *
* Department: Oden Institute and Institute for Cellular and Molecular Biology  *
* PI: Edward Marcotte                                                          *
* Project: Protein Fluorosequencing                                            *
\******************************************************************************/

// Defining symbols from header:
#include "bleach-transition.h"

// Local project headers:
#include "hmm/fit/error-model-fitter.h"
#include "hmm/step/binomial-transition.h"
#include "tensor/tensor.h"

namespace fluoroseq {

BleachTransition::BleachTransition(double q, int channel)
        : BinomialTransition(q, channel) {}

void BleachTransition::improve_fit(const Tensor& forward_tensor,
                                   const Tensor& backward_tensor,
                                   const Tensor& next_backward_tensor,
                                   int edmans,
                                   double probability,
                                   ErrorModelFitter* fitter) const {
    BinomialTransition::improve_fit(forward_tensor,
                                    backward_tensor,
                                    next_backward_tensor,
                                    edmans,
                                    probability,
                                    &fitter->p_bleach_fit);
}

}  // namespace fluoroseq
