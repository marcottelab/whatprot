/******************************************************************************\
* Author: Matthew Beauregard Smith                                             *
* Affiliation: The University of Texas at Austin                               *
* Department: Oden Institute and Institute for Cellular and Molecular Biology  *
* PI: Edward Marcotte                                                          *
* Project: Protein Fluorosequencing                                            *
\******************************************************************************/

#ifndef WHATPROT_HMM_STEP_CYCLIC_BLEACH_TRANSITION_H
#define WHATPROT_HMM_STEP_CYCLIC_BLEACH_TRANSITION_H

// Standard C++ library headers:
#include <vector>

// Local project headers:
#include "hmm/state-vector/peptide-state-vector.h"
#include "hmm/step/binomial-transition.h"
#include "parameterization/fit/sequencing-model-fitter.h"

namespace whatprot {

class CyclicBleachTransition : public BinomialTransition {
public:
    CyclicBleachTransition(double q, int channel);
    CyclicBleachTransition(const CyclicBleachTransition& other);
    virtual void improve_fit(const PeptideStateVector& forward_psv,
                             const PeptideStateVector& backward_psv,
                             const PeptideStateVector& next_backward_psv,
                             unsigned int num_edmans,
                             double probability,
                             SequencingModelFitter* fitter) const override;
};

}  // namespace whatprot

#endif  // WHATPROT_HMM_STEP_CYCLIC_BLEACH_TRANSITION_H