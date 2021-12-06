/******************************************************************************\
* Author: Matthew Beauregard Smith                                             *
* Affiliation: The University of Texas at Austin                               *
* Department: Oden Institute and Institute for Cellular and Molecular Biology  *
* PI: Edward Marcotte                                                          *
* Project: Protein Fluorosequencing                                            *
\******************************************************************************/

#ifndef WHATPROT_HMM_STEP_STUCK_DYE_STEP_H
#define WHATPROT_HMM_STEP_STUCK_DYE_STEP_H

// Local project headers:
#include "hmm/state-vector/stuck-dye-state-vector.h"
#include "hmm/step/step.h"

namespace whatprot {

class StuckDyeStep : public Step<StuckDyeStateVector> {};

}  // namespace whatprot

#endif  // WHATPROT_HMM_STEP_STUCK_DYE_STEP_H
