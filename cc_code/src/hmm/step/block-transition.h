/******************************************************************************\
* Author: Matthew Beauregard Smith                                             *
* Affiliation: The University of Texas at Austin                               *
* Department: Oden Institute and Institute for Cellular and Molecular Biology  *
* PI: Edward Marcotte                                                          *
* Project: Protein Fluorosequencing                                            *
\******************************************************************************/

#ifndef WHATPROT_HMM_STEP_BROKEN_N_TRANSITION_H
#define WHATPROT_HMM_STEP_BROKEN_N_TRANSITION_H

// Local project headers:
#include "hmm/state-vector/peptide-state-vector.h"
#include "hmm/step/peptide-step.h"
#include "parameterization/fit/parameter-fitter.h"
#include "util/kd-range.h"

namespace whatprot {

class BrokenNTransition : public PeptideStep {
public:
    BrokenNTransition(double p_block);
    virtual void prune_forward(KDRange* range, bool* allow_detached) override;
    virtual void prune_backward(KDRange* range, bool* allow_detached) override;
    virtual PeptideStateVector* forward(
            const PeptideStateVector& input,
            unsigned int* num_edmans) const override;
    virtual PeptideStateVector* backward(
            const PeptideStateVector& input,
            unsigned int* num_edmans) const override;
    void improve_fit(const PeptideStateVector& forward_psv,
                     const PeptideStateVector& backward_psv,
                     const PeptideStateVector& next_backward_psv,
                     unsigned int num_edmans,
                     double probability,
                     ParameterFitter* fitter) const;
    KDRange pruned_range;
    double p_block;
};

}  // namespace whatprot

#endif  // WHATPROT_HMM_STEP_BROKEN_N_TRANSITION_H