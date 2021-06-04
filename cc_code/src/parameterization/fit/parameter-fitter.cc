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

ParameterFitter ParameterFitter::operator+(const ParameterFitter& other) const {
    ParameterFitter result_fitter;
    result_fitter.numerator = numerator + other.numerator;
    result_fitter.denominator = denominator + other.denominator;
    return result_fitter;
}

void ParameterFitter::operator+=(const ParameterFitter& other) {
    numerator += other.numerator;
    denominator += other.denominator;
}

void ParameterFitter::operator*=(double weight_adjustment) {
    numerator *= weight_adjustment;
    denominator *= weight_adjustment;
}

}  // namespace whatprot
