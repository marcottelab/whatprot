/******************************************************************************\
* Author: Matthew Beauregard Smith                                             *
* Affiliation: The University of Texas at Austin                               *
* Department: Oden Institute and Institute for Cellular and Molecular Biology  *
* PI: Edward Marcotte                                                          *
* Project: Protein Fluorosequencing                                            *
\******************************************************************************/

// Defining symbols from header:
#include "parameter-fitter.h"

namespace whatprot {

ParameterFitter::ParameterFitter() : numerator(0.0), denominator(0.0) {}

double ParameterFitter::get() const {
    return numerator / denominator;
}

}  // namespace whatprot
