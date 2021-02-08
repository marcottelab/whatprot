/******************************************************************************\
* Author: Matthew Beauregard Smith                                             *
* Affiliation: The University of Texas at Austin                               *
* Department: Oden Institute and Institute for Cellular and Molecular Biology  *
* PI: Edward Marcotte                                                          *
* Project: Protein Fluorosequencing                                            *
\******************************************************************************/

#ifndef WHATPROT_HMM_SIMPLE_HMM_H
#define WHATPROT_HMM_SIMPLE_HMM_H

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
               const Summation& summation);

}  // namespace whatprot

#endif  // WHATPROT_HMM_SIMPLE_HMM_H
