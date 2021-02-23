/******************************************************************************\
* Author: Matthew Beauregard Smith                                             *
* Affiliation: The University of Texas at Austin                               *
* Department: Oden Institute and Institute for Cellular and Molecular Biology  *
* PI: Edward Marcotte                                                          *
* Project: Protein Fluorosequencing                                            *
\******************************************************************************/

#ifndef WHATPROT_HMM_PRECOMPUTATIONS_RADIOMETRY_PRECOMPUTATIONS_H
#define WHATPROT_HMM_PRECOMPUTATIONS_RADIOMETRY_PRECOMPUTATIONS_H

// Local project headers:
#include "common/error-model.h"
#include "common/radiometry.h"
#include "hmm/step/emission.h"

namespace whatprot {

class RadiometryPrecomputations {
public:
    RadiometryPrecomputations(const Radiometry& radiometry,
                              const ErrorModel& error_model,
                              int max_num_dyes);
    Emission emission;
};

}  // namespace whatprot

#endif  // WHATPROT_HMM_PRECOMPUTATIONS_RADIOMETRY_PRECOMPUTATIONS_H
