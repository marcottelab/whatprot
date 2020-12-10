/******************************************************************************\
* Author: Matthew Beauregard Smith                                             *
* Affiliation: The University of Texas at Austin                               *
* Department: Oden Institute and Institute for Cellular and Molecular Biology  *
* PI: Edward Marcotte                                                          *
* Project: Protein Fluorosequencing                                            *
\******************************************************************************/

#ifndef FLUOROSEQ_FWD_ALG_DETACH_TRANSITION
#define FLUOROSEQ_FWD_ALG_DETACH_TRANSITION

// Local project headers:
#include "tensor/tensor.h"

namespace fluoroseq {

class DetachTransition {
public:
    DetachTransition(double p_detach);
    void operator()(Tensor* tensor, int edmans) const;

    double p_detach;
};

}  // namespace fluoroseq

#endif  // FLUOROSEQ_FWD_ALG_DETACH_TRANSITION
