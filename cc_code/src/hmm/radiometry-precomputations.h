/******************************************************************************\
* Author: Matthew Beauregard Smith                                             *
* Affiliation: The University of Texas at Austin                               *
* Department: Oden Institute and Institute for Cellular and Molecular Biology  *
* PI: Edward Marcotte                                                          *
* Project: Protein Fluorosequencing                                            *
\******************************************************************************/

#ifndef FLUOROSEQ_HMM_RADIOMETRY_PRECOMPUTATIONS_H
#define FLUOROSEQ_HMM_RADIOMETRY_PRECOMPUTATIONS_H

// Local project headers:
#include "common/error-model.h"
#include "common/radiometry.h"
#include "hmm/emission.h"

namespace fluoroseq {

class RadiometryPrecomputations {
public:
    RadiometryPrecomputations(const Radiometry& radiometry,
                              const ErrorModel& error_model,
                              int max_num_dyes);
    Emission emission;
};

}  // namespace fluoroseq

#endif  // FLUOROSEQ_HMM_RADIOMETRY_PRECOMPUTATIONS_H
