// Author: Matthew Beauregard Smith (UT Austin)
#include "summation.h"

#include "tensor/tensor.h"

namespace fluoroseq {

double Summation::operator()(Tensor * tensor, int timestep) const {
    double sum = 0.0;
    for (int i = 0; i < timestep * tensor->strides[0]; i++) {
        sum += tensor->values[i];
    }
    return sum;
}

}  // namespace fluoroseq
