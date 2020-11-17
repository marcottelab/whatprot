/******************************************************************************\
* Author: Matthew Beauregard Smith                                             *
* Affiliation: The University of Texas at Austin                               *
* Department: Oden Institute and Institute for Cellular and Molecular Biology  *
* PI: Edward Marcotte                                                          *
* Project: Protein Fluorosequencing                                            *
\******************************************************************************/

#ifndef FLUOROSEQ_HMM_SIMPLE_HMM_H
#define FLUOROSEQ_HMM_SIMPLE_HMM_H

// Local project headers:
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
               const Summation& summation);

}  // namespace fluoroseq

#endif  // FLUOROSEQ_HMM_SIMPLE_HMM_H
