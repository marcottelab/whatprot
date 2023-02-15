/******************************************************************************\
* Author: Matthew Beauregard Smith                                             *
* Affiliation: The University of Texas at Austin                               *
* Department: Oden Institute and Institute for Cellular and Molecular Biology  *
* PI: Edward Marcotte                                                          *
* Project: Protein Fluorosequencing                                            *
\******************************************************************************/

// Defining symbols from header:
#include "cyclic-bleach-transition.h"

// Local project headers:
#include "hmm/state-vector/peptide-state-vector.h"
#include "hmm/step/binomial-transition.h"
#include "parameterization/fit/sequencing-model-fitter.h"

namespace whatprot {

CyclicBleachTransition::CyclicBleachTransition(double q, int channel)
        : BinomialTransition(q, channel) {}

CyclicBleachTransition::CyclicBleachTransition(
        const CyclicBleachTransition& other)
        : BinomialTransition(other) {}

void CyclicBleachTransition::improve_fit(
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
            &fitter->channel_fits[channel]->p_cyclic_bleach_fit);
}

}  // namespace whatprot
