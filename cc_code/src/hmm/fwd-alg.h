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
               const Summation& summation);

}  // namespace fluoroseq

#endif  // FLUOROSEQ_HMM_SIMPLE_HMM_H
