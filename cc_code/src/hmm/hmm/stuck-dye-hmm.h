/******************************************************************************\
* Author: Matthew Beauregard Smith                                             *
* Affiliation: The University of Texas at Austin                               *
* Department: Oden Institute and Institute for Cellular and Molecular Biology  *
* PI: Edward Marcotte                                                          *
* Project: Protein Fluorosequencing                                            *
\******************************************************************************/

#ifndef WHATPROT_HMM_HMM_PEPTIDE_HMM_H
#define WHATPROT_HMM_HMM_PEPTIDE_HMM_H

// Standard C++ library headers:
#include <vector>

// Local project headers:
#include "hmm/fit/error-model-fitter.h"
#include "hmm/hmm/generic-hmm.h"
#include "hmm/state-vector/stuck-dye-state-vector.h"
#include "hmm/step/step.h"

namespace whatprot {

class StuckDyeHMM : public GenericHMM<StuckDyeStateVector> {
public:
    StuckDyeHMM(int num_timesteps,
               int num_channels,
               int channel,
                const RadiometryPrecomputations& radiometry_precomputations,
                const UniversalPrecomputations& universal_precomputations);
    virtual StuckDyeStateVector create_states() const override;
};

}  // namespace whatprot

#endif  // WHATPROT_HMM_HMM_PEPTIDE_HMM_H
