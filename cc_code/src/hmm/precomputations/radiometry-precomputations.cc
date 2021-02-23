/******************************************************************************\
* Author: Matthew Beauregard Smith                                             *
* Affiliation: The University of Texas at Austin                               *
* Department: Oden Institute and Institute for Cellular and Molecular Biology  *
* PI: Edward Marcotte                                                          *
* Project: Protein Fluorosequencing                                            *
\******************************************************************************/

// Defining symbols from header:
#include "radiometry-precomputations.h"

// Local project headers:
#include "common/error-model.h"
#include "common/radiometry.h"
#include "hmm/step/emission.h"

namespace whatprot {

RadiometryPrecomputations::RadiometryPrecomputations(
        const Radiometry& radiometry,
        const ErrorModel& error_model,
        int max_num_dyes)
        : emission(radiometry, max_num_dyes, error_model.pdf()) {}

}  // namespace whatprot
