/******************************************************************************\
* Author: Matthew Beauregard Smith                                             *
* Affiliation: The University of Texas at Austin                               *
* Department: Oden Institute and Institute for Cellular and Molecular Biology  *
* PI: Edward Marcotte                                                          *
* Project: Protein Fluorosequencing                                            *
\******************************************************************************/

// Defining symbols from header:
#include "universal-precomputations.h"

// Standard C++ library headers:
#include <vector>

// Local project headers:
#include "hmm/step/bleach-transition.h"
#include "hmm/step/detach-transition.h"
#include "hmm/step/dud-transition.h"
#include "parameterization/model/sequencing-model.h"

namespace whatprot {

UniversalPrecomputations::UniversalPrecomputations(
        const SequencingModel& seq_model, unsigned int num_channels)
        : detach_transition(seq_model.p_detach),
          initial_broken_n_transition(seq_model.p_initial_break_n),
          cyclic_broken_n_transition(seq_model.p_cyclic_break_n),
          num_channels(num_channels) {
    for (unsigned int i = 0; i < num_channels; i++) {
        dud_transitions.push_back(
                new DudTransition(seq_model.channel_models[i]->p_dud, i));
        bleach_transitions.push_back(
                new BleachTransition(seq_model.channel_models[i]->p_bleach, i));
    }
}

UniversalPrecomputations::~UniversalPrecomputations() {
    for (DudTransition* step : dud_transitions) {
        delete step;
    }
    for (BleachTransition* step : bleach_transitions) {
        delete step;
    }
}

void UniversalPrecomputations::set_max_num_dyes(int max_num_dyes) {
    for (unsigned int i = 0; i < num_channels; i++) {
        dud_transitions[i]->reserve(max_num_dyes);
        bleach_transitions[i]->reserve(max_num_dyes);
    }
}

}  // namespace whatprot
