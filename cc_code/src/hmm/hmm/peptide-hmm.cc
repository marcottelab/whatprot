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
#include "hmm/step/binomial-transition.h"
#include "hmm/step/detach-transition.h"
#include "hmm/step/edman-transition.h"
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
    for (unsigned int c = 0; c < num_channels; c++) {
        steps.push_back(new DudTransition(
                universal_precomputations.dud_transitions[c]));
    }
    steps.push_back(new PeptideEmission(
            radiometry_precomputations.peptide_emissions[0]));
    for (unsigned int t = 1; t < num_timesteps; t++) {
        steps.push_back(new DetachTransition(
                universal_precomputations.detach_transition));
        for (unsigned int c = 0; c < num_channels; c++) {
            steps.push_back(new BleachTransition(
                    universal_precomputations.bleach_transitions[c]));
        }
        steps.push_back(
                new EdmanTransition(dye_seq_precomputations.edman_transition));
        steps.push_back(new PeptideEmission(
                radiometry_precomputations.peptide_emissions[t]));
    }
    tensor_shape = dye_seq_precomputations.tensor_shape;
    // Now we prune to improve efficiency when run.
    KDRange range;
    range.min = vector<unsigned int>(tensor_shape.size(), 0u);
    range.max = tensor_shape;
    bool allow_detached;
    for (unsigned int i = 0; i < steps.size(); i++) {
        steps[i]->prune_forward(&range, &allow_detached);
        if (range.is_empty()) {
            empty_range = true;
            return;
        }
    }
    for (int i = steps.size() - 1; i >= 0; i--) {
        steps[i]->prune_backward(&range, &allow_detached);
        if (range.is_empty()) {
            empty_range = true;
            return;
        }
    }
}

PeptideStateVector* PeptideHMM::create_states() const {
    return new PeptideStateVector(tensor_shape.size(), &tensor_shape[0]);
}

double PeptideHMM::probability() const {
    if (empty_range) {
        return 0.0;
    } else {
        return GenericHMM::probability();
    }
}

}  // namespace whatprot
