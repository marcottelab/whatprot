// Author: Matthew Beauregard Smith (UT Austin)
#include "detach_transition.h"

#include <algorithm>

#include "tensor/tensor.h"

namespace fluoroseq {

namespace {
using std::max;
}  // namespace

DetachTransition::DetachTransition(double p_detach, int max_failed_edmans)
        : p_detach(p_detach),
          max_failed_edmans(max_failed_edmans) {}

void DetachTransition::operator()(Tensor* tensor, int edmans) const {
    int i_min = max(0, edmans - max_failed_edmans) * tensor->strides[0];
    int i_max = (edmans + 1) * tensor->strides[0];
    double sum = 0.0;
    for (int i = i_min; i < i_max; i++) {
        double value = tensor->values[i];
        tensor->values[i] = value * (1 - p_detach);
        sum += value;
    }
    tensor->values[edmans * tensor->strides[0]] += p_detach * sum;
}

}  // namespace fluoroseq
