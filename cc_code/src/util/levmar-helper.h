/******************************************************************************\
* Author: Matthew Beauregard Smith                                             *
* Affiliation: Erisyon, Inc.                                                   *
* Project: Protein Fluorosequencing                                            *
\******************************************************************************/

#ifndef WHATPROT_UTIL_LEVMAR_HELPER_H
#define WHATPROT_UTIL_LEVMAR_HELPER_H

// Standard C++ library headers:
#include <vector>

namespace whatprot {

void offset_exponential(double* p, double* y, int m, int n, void* data);

void jacobian_of_offset_exponential(
        double* p, double* jac, int m, int n, void* data);

void least_squares_fit_of_offset_exponential(const std::vector<double>& y,
                                             const std::vector<double>& w,
                                             const std::vector<bool>& hold,
                                             std::vector<double>* p);

}  // namespace whatprot

#endif  // WHATPROT_UTIL_LEVMAR_HELPER_H
