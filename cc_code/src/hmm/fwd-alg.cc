/******************************************************************************\
* Author: Matthew Beauregard Smith                                             *
* Affiliation: The University of Texas at Austin                               *
* Department: Oden Institute and Institute for Cellular and Molecular Biology  *
* PI: Edward Marcotte                                                          *
* Project: Protein Fluorosequencing                                            *
\******************************************************************************/

// Defining symbols from header:
#include "fwd-alg.h"

// Local project headers:
#include "hmm/binomial-transition.h"
#include "hmm/detach-transition.h"
#include "hmm/edman-transition.h"
#include "hmm/emission.h"
#include "hmm/initialization.h"
#include "hmm/summation.h"
#include "tensor/tensor.h"

namespace fluoroseq {

double fwd_alg(Tensor* states,
               int num_timesteps,
               int num_channels,
               const Initialization& initialization,
               const Emission& emission,
               const DetachTransition& detach_transition,
               const BinomialTransition& dud_transition,
               const BinomialTransition& bleach_transition,
               const EdmanTransition& edman_transition,
               const Summation& summation) {
    initialization.forward(states);
    for (int c = 0; c < num_channels; c++) {
        dud_transition.forward(*states, c, 0, states);
    }
    emission.forward(*states, 0, states);
    for (int t = 1; t < num_timesteps; t++) {
        detach_transition.forward(*states, t - 1, states);
        for (int c = 0; c < num_channels; c++) {
            bleach_transition.forward(*states, c, t - 1, states);
        }
        edman_transition.forward(*states, t - 1, states);
        emission.forward(*states, t, states);
    }
    return summation(*states, num_timesteps - 1);
}

}  // namespace fluoroseq
