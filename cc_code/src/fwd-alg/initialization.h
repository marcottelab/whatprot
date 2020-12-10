/******************************************************************************\
* Author: Matthew Beauregard Smith                                             *
* Affiliation: The University of Texas at Austin                               *
* Department: Oden Institute and Institute for Cellular and Molecular Biology  *
* PI: Edward Marcotte                                                          *
* Project: Protein Fluorosequencing                                            *
\******************************************************************************/

#ifndef FLUOROSEQ_FWD_ALG_INITIALIZATION
#define FLUOROSEQ_FWD_ALG_INITIALIZATION

// Local project headers:
#include "tensor/tensor.h"

namespace fluoroseq {

class Initialization {
public:
    void operator()(Tensor* tensor) const;
};

}  // namespace fluoroseq

#endif  // FLUOROSEQ_FWD_ALG_INITIALIZATION