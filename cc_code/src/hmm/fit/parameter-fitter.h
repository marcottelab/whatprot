/******************************************************************************\
* Author: Matthew Beauregard Smith                                             *
* Affiliation: The University of Texas at Austin                               *
* Department: Oden Institute and Institute for Cellular and Molecular Biology  *
* PI: Edward Marcotte                                                          *
* Project: Protein Fluorosequencing                                            *
\******************************************************************************/

#ifndef WHATPROT_HMM_FIT_PARAMETER_FITTER_H
#define WHATPROT_HMM_FIT_PARAMETER_FITTER_H

namespace whatprot {

class ParameterFitter {
public:
    ParameterFitter();
    double get() const;
    double numerator;
    double denominator;
};

}  // namespace whatprot

#endif  // WHATPROT_HMM_FIT_PARAMETER_FITTER_H
