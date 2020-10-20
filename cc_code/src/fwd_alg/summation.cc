// Author: Matthew Beauregard Smith (UT Austin)
#include "summation.h"

#include "tensor/tensor.h"

namespace fluoroseq {

Summation::Summation() {}

double Summation::operator()(Tensor * tensor, int timestep) const {
    double sum = 0.0;
    for (int i = 0; i < (timestep + 1) * tensor->strides[0]; i++) {
        sum += tensor->values[i];
    }
    return sum;
}

}  // namespace fluoroseq
