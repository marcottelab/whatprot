/******************************************************************************\
* Author: Matthew Beauregard Smith                                             *
* Affiliation: The University of Texas at Austin                               *
* Department: Oden Institute and Institute for Cellular and Molecular Biology  *
* PI: Edward Marcotte                                                          *
* Project: Protein Fluorosequencing                                            *
\******************************************************************************/

#ifndef WHATPROT_HMM_STEP_START_H
#define WHATPROT_HMM_STEP_START_H

// Local project headers:
#include "hmm/fit/error-model-fitter.h"
#include "hmm/state-vector/peptide-state-vector.h"
#include "hmm/step/step.h"

namespace whatprot {

class Start : public Step<PeptideStateVector> {
public:
    virtual void forward(int* num_edmans, PeptideStateVector* psv) const override;
    virtual void backward(const PeptideStateVector& input, int* num_edmans,
                          PeptideStateVector* output) const override;
    virtual void improve_fit(const PeptideStateVector& forward_psv,
                             const PeptideStateVector& backward_psv,
                             const PeptideStateVector& next_backward_psv, int num_edmans,
                             double probability,
                             ErrorModelFitter* fitter) const override;
};

}  // namespace whatprot

#endif  // WHATPROT_HMM_STEP_START_H