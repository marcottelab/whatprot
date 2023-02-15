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
#include "hmm/step/cyclic-bleach-transition.h"
#include "hmm/step/cyclic-broken-n-transition.h"
#include "hmm/step/cyclic-detach-transition.h"
#include "hmm/step/dud-transition.h"
#include "hmm/step/initial-bleach-transition.h"
#include "hmm/step/initial-broken-n-transition.h"
#include "hmm/step/initial-detach-transition.h"
#include "parameterization/model/sequencing-model.h"

namespace whatprot {

class UniversalPrecomputations {
public:
    UniversalPrecomputations(const SequencingModel& seq_model,
                             unsigned int num_channels);
    ~UniversalPrecomputations();
    void set_max_num_dyes(int max_num_dyes);
    InitialDetachTransition initial_detach_transition;
    CyclicDetachTransition cyclic_detach_transition;
    InitialBrokenNTransition initial_broken_n_transition;
    CyclicBrokenNTransition cyclic_broken_n_transition;
    std::vector<DudTransition*> dud_transitions;
    std::vector<InitialBleachTransition*> initial_bleach_transitions;
    std::vector<CyclicBleachTransition*> cyclic_bleach_transitions;
    unsigned int num_channels;
};

}  // namespace whatprot

#endif  // WHATPROT_HMM_PRECOMPUTATIONS_UNIVERSAL_PRECOMPUTATIONS_H
