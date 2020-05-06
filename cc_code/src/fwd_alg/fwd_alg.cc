// Author: Matthew Beauregard Smith
#include "fwd_alg.h"

#include "fwd_alg/binomial_transition.h"
#include "fwd_alg/detach_transition.h"
#include "fwd_alg/edman_transition.h"
#include "fwd_alg/emission.h"
#include "fwd_alg/initialization.h"
#include "fwd_alg/summation.h"
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
        edman_transition(states, t);
        emission(states, t);
    }
    return summation(states, num_timesteps);
}

}  // namespace fluoroseq
