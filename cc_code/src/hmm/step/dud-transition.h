/******************************************************************************\
* Author: Matthew Beauregard Smith                                             *
* Affiliation: The University of Texas at Austin                               *
* Department: Oden Institute and Institute for Cellular and Molecular Biology  *
* PI: Edward Marcotte                                                          *
* Project: Protein Fluorosequencing                                            *
\******************************************************************************/

#ifndef FLUOROSEQ_HMM_STEP_DUD_TRANSITION_H
#define FLUOROSEQ_HMM_STEP_DUD_TRANSITION_H

// Standard C++ library headers:
#include <vector>

// Local project headers:
#include "hmm/fit/error-model-fitter.h"
#include "hmm/step/binomial-transition.h"
#include "tensor/tensor.h"

namespace fluoroseq {

class DudTransition : public BinomialTransition {
public:
    DudTransition(double q, int channel);
    virtual void improve_fit(const Tensor& forward_tensor,
                             const Tensor& backward_tensor,
                             const Tensor& next_backward_tensor,
                             int edmans,
                             double probability,
                             ErrorModelFitter* fitter) const override;
};

}  // namespace fluoroseq

#endif  // FLUOROSEQ_HMM_STEP_DUD_TRANSITION_H