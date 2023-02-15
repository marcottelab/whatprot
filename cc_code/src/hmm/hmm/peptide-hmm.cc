/******************************************************************************\
* Author: Matthew Beauregard Smith                                             *
* Affiliation: The University of Texas at Austin                               *
* Department: Oden Institute and Institute for Cellular and Molecular Biology  *
* PI: Edward Marcotte                                                          *
* Project: Protein Fluorosequencing                                            *
\******************************************************************************/

// Defining symbols from header:
#include "peptide-hmm.h"

// Standard C++ library headers:
#include <vector>

// Local project headers:
#include "hmm/precomputations/dye-seq-precomputations.h"
#include "hmm/precomputations/radiometry-precomputations.h"
#include "hmm/precomputations/universal-precomputations.h"
#include "hmm/state-vector/peptide-state-vector.h"
#include "hmm/step/cyclic-bleach-transition.h"
#include "hmm/step/cyclic-broken-n-transition.h"
#include "hmm/step/cyclic-detach-transition.h"
#include "hmm/step/edman-transition.h"
#include "hmm/step/initial-bleach-transition.h"
#include "hmm/step/initial-broken-n-transition.h"
#include "hmm/step/initial-detach-transition.h"
#include "hmm/step/peptide-emission.h"

namespace whatprot {

namespace {
using std::vector;
}

PeptideHMM::PeptideHMM(
        unsigned int num_timesteps,
        unsigned int num_channels,
        const DyeSeqPrecomputations& dye_seq_precomputations,
        const RadiometryPrecomputations& radiometry_precomputations,
        const UniversalPrecomputations& universal_precomputations)
        : GenericHMM(num_timesteps), empty_range(false) {
    steps.push_back(new InitialBrokenNTransition(
            universal_precomputations.initial_broken_n_transition));
    for (unsigned int c = 0; c < num_channels; c++) {
        steps.push_back(new DudTransition(
                *universal_precomputations.dud_transitions[c]));
    }
    steps.push_back(new PeptideEmission(
            *radiometry_precomputations.peptide_emissions[0]));
    if (num_timesteps > 1) {
        steps.push_back(new InitialDetachTransition(
                universal_precomputations.initial_detach_transition));
        for (unsigned int c = 0; c < num_channels; c++) {
            steps.push_back(new InitialBleachTransition(
                    *universal_precomputations.initial_bleach_transitions[c]));
        }
    }
    for (unsigned int t = 1; t < num_timesteps; t++) {
        steps.push_back(new CyclicBrokenNTransition(
                universal_precomputations.cyclic_broken_n_transition));
        steps.push_back(new CyclicDetachTransition(
                universal_precomputations.cyclic_detach_transition));
        for (unsigned int c = 0; c < num_channels; c++) {
            steps.push_back(new CyclicBleachTransition(
                    *universal_precomputations.cyclic_bleach_transitions[c]));
        }
        steps.push_back(
                new EdmanTransition(dye_seq_precomputations.edman_transition));
        steps.push_back(new PeptideEmission(
                *radiometry_precomputations.peptide_emissions[t]));
    }
    // Now we prune to improve efficiency when run.
    KDRange range;
    range.min = vector<unsigned int>(
            dye_seq_precomputations.tensor_shape.size(), 0u);
    range.max = dye_seq_precomputations.tensor_shape;
    bool allow_detached;
    for (unsigned int i = 0; i < steps.size(); i++) {
        steps[i]->prune_forward(&range, &allow_detached);
        if (range.is_empty()) {
            empty_range = true;
            return;
        }
    }
    backward_range = range;
    for (int i = steps.size() - 1; i >= 0; i--) {
        steps[i]->prune_backward(&range, &allow_detached);
        if (range.is_empty()) {
            empty_range = true;
            return;
        }
    }
    forward_range = range;
}

PeptideStateVector* PeptideHMM::create_states_forward() const {
    return new PeptideStateVector(forward_range);
}

PeptideStateVector* PeptideHMM::create_states_backward() const {
    return new PeptideStateVector(backward_range);
}

double PeptideHMM::probability() const {
    if (empty_range) {
        return 0.0;
    } else {
        return GenericHMM::probability();
    }
}

}  // namespace whatprot
