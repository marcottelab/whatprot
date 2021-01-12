/******************************************************************************\
* Author: Matthew Beauregard Smith                                             *
* Affiliation: The University of Texas at Austin                               *
* Department: Oden Institute and Institute for Cellular and Molecular Biology  *
* PI: Edward Marcotte                                                          *
* Project: Protein Fluorosequencing                                            *
\******************************************************************************/

#ifndef FLUOROSEQ_HMM_DETACH_TRANSITION
#define FLUOROSEQ_HMM_DETACH_TRANSITION

// Local project headers:
#include "tensor/tensor.h"

namespace fluoroseq {

class DetachTransition {
public:
    DetachTransition(double p_detach);
    void forward(const Tensor& input, int edmans, Tensor* output) const;

    double p_detach;
};

}  // namespace fluoroseq

#endif  // FLUOROSEQ_HMM_DETACH_TRANSITION
