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
#include "fwd-alg/binomial-transition.h"
#include "fwd-alg/detach-transition.h"
#include "fwd-alg/edman-transition.h"
#include "fwd-alg/emission.h"
#include "fwd-alg/initialization.h"
#include "fwd-alg/summation.h"
#include "tensor/tensor.h"

namespace whatprot {

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
    initialization(states);
    for (int c = 0; c < num_channels; c++) {
        dud_transition(states, c, 0);
    }
    emission(states, 0);
    for (int t = 1; t < num_timesteps; t++) {
        detach_transition(states, t - 1);
        for (int c = 0; c < num_channels; c++) {
            bleach_transition(states, c, t - 1);
        }
        edman_transition(states, t - 1);
        emission(states, t);
    }
    return summation(states, num_timesteps - 1);
}

}  // namespace whatprot
