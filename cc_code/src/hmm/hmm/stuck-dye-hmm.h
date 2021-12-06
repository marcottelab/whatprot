/******************************************************************************\
* Author: Matthew Beauregard Smith                                             *
* Affiliation: The University of Texas at Austin                               *
* Department: Oden Institute and Institute for Cellular and Molecular Biology  *
* PI: Edward Marcotte                                                          *
* Project: Protein Fluorosequencing                                            *
\******************************************************************************/

#ifndef WHATPROT_HMM_HMM_STUCK_DYE_HMM_H
#define WHATPROT_HMM_HMM_STUCK_DYE_HMM_H

// Standard C++ library headers:
#include <vector>

// Local project headers:
#include "hmm/hmm/generic-hmm.h"
#include "hmm/precomputations/radiometry-precomputations.h"
#include "hmm/precomputations/universal-precomputations.h"
#include "hmm/state-vector/stuck-dye-state-vector.h"
#include "hmm/step/stuck-dye-step.h"
#include "parameterization/fit/sequencing-model-fitter.h"

namespace whatprot {

class StuckDyeHMM : public GenericHMM<StuckDyeStateVector, StuckDyeStep> {
public:
    StuckDyeHMM(unsigned int num_timesteps,
                unsigned int num_channels,
                unsigned int channel,
                const RadiometryPrecomputations& radiometry_precomputations,
                const UniversalPrecomputations& universal_precomputations);
    virtual StuckDyeStateVector create_states() const override;
};

}  // namespace whatprot

#endif  // WHATPROT_HMM_HMM_STUCK_DYE_HMM_H
