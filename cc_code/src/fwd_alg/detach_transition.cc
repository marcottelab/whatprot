// Author: Matthew Beauregard Smith (UT Austin)
#include "detach_transition.h"

#include "tensor/tensor.h"

namespace fluoroseq {

DetachTransition::DetachTransition(double p_detach, int max_failed_edmans)
        : p_detach(p_detach),
          max_failed_edmans(max_failed_edmans) {}

void DetachTransition::operator()(Tensor* tensor, int edmans) const {
    int e_min = edmans - max_failed_edmans;
    if (e_min < 0) {
        e_min = 0;
    }
    for (int e = e_min; e < edmans + 1; e++) {
        double sum = 0.0;
        for (int i = 1; i < tensor->strides[0]; i++) {
            double value = tensor->values[e * tensor->strides[0] + i];
            tensor->values[e * tensor->strides[0] + i] = value * (1 - p_detach);
            sum += value;
        }
        tensor->values[e * tensor->strides[0]] += p_detach * sum;
    }
}

}  // namespace fluoroseq
