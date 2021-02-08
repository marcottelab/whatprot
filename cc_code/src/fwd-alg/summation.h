/******************************************************************************\
* Author: Matthew Beauregard Smith                                             *
* Affiliation: The University of Texas at Austin                               *
* Department: Oden Institute and Institute for Cellular and Molecular Biology  *
* PI: Edward Marcotte                                                          *
* Project: Protein Fluorosequencing                                            *
\******************************************************************************/

#ifndef WHATPROT_FWD_ALG_SUMMATION_H
#define WHATPROT_FWD_ALG_SUMMATION_H

// Local project headers:
#include "tensor/tensor.h"

namespace whatprot {

class Summation {
public:
    Summation();
    double operator()(Tensor* tensor, int timestep) const;
};

}  // namespace whatprot

#endif  // WHATPROT_FWD_ALG_SUMMATION_H