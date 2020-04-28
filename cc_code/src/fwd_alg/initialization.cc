// Author: Matthew Beauregard Smith (UT Austin)
#include "initialization.h"

#include "tensor/tensor.h"

namespace fluoroseq {

void Initialization::operator()(Tensor* tensor) const {
    for (int i = 0; i < tensor->strides[0] - 1; i++) {
        tensor->values[i] = 0.0;
    }
    tensor->values[tensor->strides[0] - 1] = 1.0;
}

}  // namespace fluoroseq
