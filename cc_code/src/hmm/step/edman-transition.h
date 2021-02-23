/******************************************************************************\
* Author: Matthew Beauregard Smith                                             *
* Affiliation: The University of Texas at Austin                               *
* Department: Oden Institute and Institute for Cellular and Molecular Biology  *
* PI: Edward Marcotte                                                          *
* Project: Protein Fluorosequencing                                            *
\******************************************************************************/

#ifndef WHATPROT_HMM_STEP_EDMAN_TRANSITION_H
#define WHATPROT_HMM_STEP_EDMAN_TRANSITION_H

// Local project headers:
#include "common/dye-track.h"
#include "hmm/fit/error-model-fitter.h"
#include "hmm/step/step.h"
#include "tensor/tensor.h"

namespace whatprot {

class EdmanTransition : public Step {
public:
    EdmanTransition(double p_edman_failure,
                    const DyeSeq& dye_seq,
                    const DyeTrack& dye_track);
    virtual void forward(int* edmans,
                         Tensor* tsr) const override;
    virtual void backward(const Tensor& input,
                          int* edmans,
                          Tensor* output) const override;
    virtual void improve_fit(const Tensor& forward_tensor,
                             const Tensor& backward_tensor,
                             const Tensor& next_backward_tensor,
                             int edmans,
                             double probability,
                             ErrorModelFitter* fitter) const override;

    DyeSeq dye_seq;
    DyeTrack dye_track;
    double p_edman_failure;
};

}  // namespace whatprot

#endif  // WHATPROT_HMM_STEP_EDMAN_TRANSITION_H