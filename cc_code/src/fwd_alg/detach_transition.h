// Author: Matthew Beauregard Smith (UT Austin)
#ifndef FLUOROSEQ_FWD_ALG_DETACH_TRANSITION
#define FLUOROSEQ_FWD_ALG_DETACH_TRANSITION

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
