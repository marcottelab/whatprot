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

namespace whatprot {

enum DistributionType {
    // NORMAL,
    LOGNORMAL,
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
