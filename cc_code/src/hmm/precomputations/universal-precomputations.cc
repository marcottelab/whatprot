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
#include "hmm/step/cyclic-bleach-transition.h"
#include "hmm/step/cyclic-broken-n-transition.h"
#include "hmm/step/cyclic-detach-transition.h"
#include "hmm/step/dud-transition.h"
#include "hmm/step/initial-bleach-transition.h"
#include "hmm/step/initial-broken-n-transition.h"
#include "hmm/step/initial-detach-transition.h"
#include "parameterization/model/sequencing-model.h"

namespace whatprot {

UniversalPrecomputations::UniversalPrecomputations(
        const SequencingModel& seq_model, unsigned int num_channels)
        : initial_detach_transition(seq_model.p_initial_detach),
          cyclic_detach_transition(seq_model.p_cyclic_detach),
          initial_broken_n_transition(seq_model.p_initial_break_n),
          cyclic_broken_n_transition(seq_model.p_cyclic_break_n),
          num_channels(num_channels) {
    for (unsigned int i = 0; i < num_channels; i++) {
        dud_transitions.push_back(
                new DudTransition(seq_model.channel_models[i]->p_dud, i));
        initial_bleach_transitions.push_back(new InitialBleachTransition(
                seq_model.channel_models[i]->p_initial_bleach, i));
        cyclic_bleach_transitions.push_back(new CyclicBleachTransition(
                seq_model.channel_models[i]->p_cyclic_bleach, i));
    }
}

UniversalPrecomputations::~UniversalPrecomputations() {
    for (DudTransition* step : dud_transitions) {
        delete step;
    }
    for (InitialBleachTransition* step : initial_bleach_transitions) {
        delete step;
    }
    for (CyclicBleachTransition* step : cyclic_bleach_transitions) {
        delete step;
    }
}

void UniversalPrecomputations::set_max_num_dyes(int max_num_dyes) {
    for (unsigned int i = 0; i < num_channels; i++) {
        dud_transitions[i]->reserve(max_num_dyes);
        initial_bleach_transitions[i]->reserve(max_num_dyes);
        cyclic_bleach_transitions[i]->reserve(max_num_dyes);
    }
}

}  // namespace whatprot
