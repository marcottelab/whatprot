/******************************************************************************\
* Author: Matthew Beauregard Smith                                             *
* Affiliation: The University of Texas at Austin                               *
* Department: Oden Institute and Institute for Cellular and Molecular Biology  *
* PI: Edward Marcotte                                                          *
* Project: Protein Fluorosequencing                                            *
\******************************************************************************/

#ifndef WHATPROT_HMM_PRECOMPUTATIONS_UNIVERSAL_PRECOMPUTATIONS_H
#define WHATPROT_HMM_PRECOMPUTATIONS_UNIVERSAL_PRECOMPUTATIONS_H

// Standard C++ library headers:
#include <vector>

// Local project headers:
#include "hmm/step/bleach-transition.h"
#include "hmm/step/cyclic-block-transition.h"
#include "hmm/step/detach-transition.h"
#include "hmm/step/dud-transition.h"
#include "hmm/step/initial-block-transition.h"
#include "parameterization/model/sequencing-model.h"

namespace whatprot {

class UniversalPrecomputations {
public:
    UniversalPrecomputations(const SequencingModel& seq_model,
                             unsigned int num_timesteps,
                             unsigned int num_channels);
    ~UniversalPrecomputations();
    void set_max_num_dyes(int max_num_dyes);
    std::vector<DetachTransition*> detach_transitions;
    InitialBrokenNTransition initial_broken_n_transition;
    CyclicBrokenNTransition cyclic_broken_n_transition;
    std::vector<DudTransition*> dud_transitions;
    std::vector<BleachTransition*> bleach_transitions;
    unsigned int num_channels;
};

}  // namespace whatprot

#endif  // WHATPROT_HMM_PRECOMPUTATIONS_UNIVERSAL_PRECOMPUTATIONS_H
