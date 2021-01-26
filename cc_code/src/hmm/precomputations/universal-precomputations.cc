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
#include "common/error-model.h"
#include "hmm/step/binomial-transition.h"
#include "hmm/step/detach-transition.h"

namespace fluoroseq {

UniversalPrecomputations::UniversalPrecomputations(
        const ErrorModel& error_model, int num_channels)
        : detach_transition(error_model.p_detach), num_channels(num_channels) {
    for (int i = 0; i < num_channels; i++) {
        dud_transitions.emplace_back(error_model.p_dud, i);
        bleach_transitions.emplace_back(error_model.p_bleach, i);
    }
}

void UniversalPrecomputations::set_max_num_dyes(int max_num_dyes) {
    for (int i = 0; i < num_channels; i++) {
        dud_transitions[i].reserve(max_num_dyes);
        bleach_transitions[i].reserve(max_num_dyes);
    }
}

}  // namespace fluoroseq
