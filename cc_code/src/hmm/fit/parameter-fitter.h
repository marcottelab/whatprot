/******************************************************************************\
* Author: Matthew Beauregard Smith                                             *
* Affiliation: The University of Texas at Austin                               *
* Department: Oden Institute and Institute for Cellular and Molecular Biology  *
* PI: Edward Marcotte                                                          *
* Project: Protein Fluorosequencing                                            *
\******************************************************************************/

#ifndef FLUOROSEQ_HMM_FIT_PARAMETER_FITTER_H
#define FLUOROSEQ_HMM_FIT_PARAMETER_FITTER_H

namespace fluoroseq {

class ParameterFitter {
public:
    ParameterFitter();
    double get() const;
    double numerator;
    double denominator;
};

}  // namespace fluoroseq

#endif  // FLUOROSEQ_HMM_FIT_PARAMETER_FITTER_H
