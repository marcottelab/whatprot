// Author: Matthew Beauregard Smith (UT Austin)
#include "summation.h"

#include "tensor/tensor.h"

namespace fluoroseq {

Summation::Summation(int max_edman_failures)
        : max_edman_failures(max_edman_failures) {}

double Summation::operator()(Tensor * tensor, int timestep) const {
    double sum = 0.0;
    int i_min = timestep - max_edman_failures;
    if (i_min < 0) {
        i_min = 0;
    }
    for (int i = i_min; i < timestep * tensor->strides[0]; i++) {
        sum += tensor->values[i];
    }
    return sum;
}

}  // namespace fluoroseq
