/******************************************************************************\
* Author: Matthew Beauregard Smith                                             *
* Affiliation: The University of Texas at Austin                               *
* Department: Oden Institute and Institute for Cellular and Molecular Biology  *
* PI: Edward Marcotte                                                          *
* Project: Protein Fluorosequencing                                            *
\******************************************************************************/

#ifndef WHATPROT_COMMON_ERROR_MODEL_H
#define WHATPROT_COMMON_ERROR_MODEL_H

// Standard C++ library headers:
#include <functional>
#include <string>

namespace whatprot {

enum DistributionType {
    NORMAL,
    LOGNORMAL,
    OVERRIDE,  // Intended for testing. Always returns 1.0 from distribution.
};

class ErrorModel {
public:
    ErrorModel(double p_edman_failure,
               double p_detach,
               double p_bleach,
               double p_dud,
               DistributionType distribution_type,
               double mu,
               double sigma);
    std::function<double(double, int)> pdf() const;
    double relative_distance(const ErrorModel& error_model) const;
    std::string debug_string() const;

    double p_edman_failure;
    double p_detach;
    double p_bleach;
    double p_dud;
    DistributionType distribution_type;
    double mu;
    double sigma;
};

}  // namespace whatprot

#endif  // WHATPROT_COMMON_ERROR_MODEL_H
