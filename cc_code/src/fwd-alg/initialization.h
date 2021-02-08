/******************************************************************************\
* Author: Matthew Beauregard Smith                                             *
* Affiliation: The University of Texas at Austin                               *
* Department: Oden Institute and Institute for Cellular and Molecular Biology  *
* PI: Edward Marcotte                                                          *
* Project: Protein Fluorosequencing                                            *
\******************************************************************************/

#ifndef WHATPROT_FWD_ALG_INITIALIZATION
#define WHATPROT_FWD_ALG_INITIALIZATION

// Local project headers:
#include "tensor/tensor.h"

namespace whatprot {

class Initialization {
public:
    void operator()(Tensor* tensor) const;
};

}  // namespace whatprot

#endif  // WHATPROT_FWD_ALG_INITIALIZATION