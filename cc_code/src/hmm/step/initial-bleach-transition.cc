/******************************************************************************\
* Author: Matthew Beauregard Smith                                             *
* Affiliation: The University of Texas at Austin                               *
* Department: Oden Institute and Institute for Cellular and Molecular Biology  *
* PI: Edward Marcotte                                                          *
* Project: Protein Fluorosequencing                                            *
\******************************************************************************/

// Defining symbols from header:
#include "initial-bleach-transition.h"

// Local project headers:
#include "hmm/state-vector/peptide-state-vector.h"
#include "hmm/step/binomial-transition.h"
#include "parameterization/fit/sequencing-model-fitter.h"

namespace whatprot {

InitialBleachTransition::InitialBleachTransition(double q, int channel)
        : BinomialTransition(q, channel) {}

InitialBleachTransition::InitialBleachTransition(
        const InitialBleachTransition& other)
        : BinomialTransition(other) {}

void InitialBleachTransition::improve_fit(
        const PeptideStateVector& forward_psv,
        const PeptideStateVector& backward_psv,
        const PeptideStateVector& next_backward_psv,
        unsigned int num_edmans,
        double probability,
        SequencingModelFitter* fitter) const {
    BinomialTransition::improve_fit(
            forward_psv,
            backward_psv,
            next_backward_psv,
            num_edmans,
            probability,
            &fitter->channel_fits[channel]->p_initial_bleach_fit);
}

}  // namespace whatprot
