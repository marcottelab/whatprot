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
#include "hmm/fit/error-model-fitter.h"
#include "hmm/hmm/generic-hmm.h"
#include "hmm/precomputations/radiometry-precomputations.h"
#include "hmm/precomputations/universal-precomputations.h"
#include "hmm/state-vector/stuck-dye-state-vector.h"
#include "hmm/step/step.h"

namespace whatprot {

StuckDyeHMM::StuckDyeHMM(
        int num_timesteps,
        int num_channels,
        int channel,
        const RadiometryPrecomputations& radiometry_precomputations,
        const UniversalPrecomputations& universal_precomputations)
        : GenericHMM(num_timesteps) {
    steps.push_back(&radiometry_precomputations.stuck_dye_emissions[channel]);
    for (int i = 1; i < num_timesteps; i++) {
        steps.push_back(&universal_precomputations.stuck_dye_transition);
        steps.push_back(
                &radiometry_precomputations.stuck_dye_emissions[channel]);
    }
}

StuckDyeStateVector StuckDyeHMM::create_states() const {
    return StuckDyeStateVector();
}

}  // namespace whatprot
