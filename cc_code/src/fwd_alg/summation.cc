// Author: Matthew Beauregard Smith (UT Austin)
#include "summation.h"

#include <algorithm>

#include "tensor/tensor.h"

namespace fluoroseq {

namespace {
using std::max;
}  // namespace

Summation::Summation(int max_edman_failures)
        : max_edman_failures(max_edman_failures) {}

double Summation::operator()(Tensor * tensor, int timestep) const {
    double sum = 0.0;
    int i_min = max(0, timestep - max_edman_failures);
    for (int i = i_min; i < (timestep + 1) * tensor->strides[0]; i++) {
        sum += tensor->values[i];
    }
    return sum;
}

}  // namespace fluoroseq
