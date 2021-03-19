/******************************************************************************\
* Author: Matthew Beauregard Smith                                             *
* Affiliation: The University of Texas at Austin                               *
* Department: Oden Institute and Institute for Cellular and Molecular Biology  *
* PI: Edward Marcotte                                                          *
* Project: Protein Fluorosequencing                                            *
\******************************************************************************/

// Defining symbols from header:
#include "stuck-dye-state-vector.h"

// Standard C++ library headers:
#include <vector>

namespace whatprot {

StuckDyeStateVector::StuckDyeStateVector(int order, const int* shape)
        : dye(0), no_dye(0) {}

double StuckDyeStateVector::sum() const {
    return dye + no_dye;
}

double StuckDyeStateVector::source() const {
    return dye;
}

}  // namespace whatprot
