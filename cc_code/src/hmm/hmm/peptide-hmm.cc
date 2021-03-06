/******************************************************************************\
* Author: Matthew Beauregard Smith                                             *
* Affiliation: The University of Texas at Austin                               *
* Department: Oden Institute and Institute for Cellular and Molecular Biology  *
* PI: Edward Marcotte                                                          *
* Project: Protein Fluorosequencing                                            *
\******************************************************************************/

// Defining symbols from header:
#include "peptide-hmm.h"

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

PeptideHMM::PeptideHMM(
        int num_timesteps,
        int num_channels,
        const DyeSeqPrecomputations& dye_seq_precomputations,
        const RadiometryPrecomputations& radiometry_precomputations,
        const UniversalPrecomputations& universal_precomputations)
        : GenericHMM(num_timesteps) {
    for (int c = 0; c < num_channels; c++) {
        steps.push_back(&universal_precomputations.dud_transitions[c]);
    }
    steps.push_back(&radiometry_precomputations.peptide_emission);
    for (int t = 1; t < num_timesteps; t++) {
        steps.push_back(&universal_precomputations.detach_transition);
        for (int c = 0; c < num_channels; c++) {
            steps.push_back(&universal_precomputations.bleach_transitions[c]);
        }
        steps.push_back(&dye_seq_precomputations.edman_transition);
        steps.push_back(&radiometry_precomputations.peptide_emission);
    }
    tensor_shape = dye_seq_precomputations.tensor_shape;
}

PeptideStateVector PeptideHMM::create_states() const {
    return PeptideStateVector(tensor_shape.size(), &tensor_shape[0]);
}

}  // namespace whatprot
