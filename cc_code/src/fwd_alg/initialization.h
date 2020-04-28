// Author: Matthew Beauregard Smith (UT Austin)
#ifndef FLUOROSEQ_FWD_ALG_INITIALIZATION
#define FLUOROSEQ_FWD_ALG_INITIALIZATION

#include "tensor/tensor.h"

namespace fluoroseq {

class Initialization {
public:
    void operator()(Tensor* tensor) const;
};

}  // namespace fluoroseq

#endif  // FLUOROSEQ_FWD_ALG_INITIALIZATION