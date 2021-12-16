/******************************************************************************\
* Author: Matthew Beauregard Smith                                             *
* Affiliation: The University of Texas at Austin                               *
* Department: Oden Institute and Institute for Cellular and Molecular Biology  *
* PI: Edward Marcotte                                                          *
* Project: Protein Fluorosequencing                                            *
\******************************************************************************/

// Defining symbols from header:
#include "stuck-dye-hmm.h"

// Standard C++ library headers:
#include <vector>

// Local project headers:
#include "hmm/hmm/generic-hmm.h"
#include "hmm/precomputations/radiometry-precomputations.h"
#include "hmm/precomputations/universal-precomputations.h"
#include "hmm/state-vector/stuck-dye-state-vector.h"
#include "hmm/step/step.h"
#include "parameterization/fit/sequencing-model-fitter.h"

namespace whatprot {

StuckDyeHMM::StuckDyeHMM(
        unsigned int num_timesteps,
        unsigned int num_channels,
        unsigned int channel,
        const RadiometryPrecomputations& radiometry_precomputations,
        const UniversalPrecomputations& universal_precomputations)
        : GenericHMM(num_timesteps) {
    steps.push_back(new StuckDyeEmission(
            radiometry_precomputations.stuck_dye_emissions[channel]));
    for (unsigned int i = 1; i < num_timesteps; i++) {
        steps.push_back(new StuckDyeTransition(
                universal_precomputations.stuck_dye_transitions[channel]));
        steps.push_back(new StuckDyeEmission(
                radiometry_precomputations.stuck_dye_emissions[channel]));
    }
}

StuckDyeStateVector* StuckDyeHMM::create_states() const {
    return new StuckDyeStateVector();
}

}  // namespace whatprot
