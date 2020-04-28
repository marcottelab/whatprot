// Author: Matthew Beauregard Smith (UT Austin)
#include "detach_transition.h"

#include "tensor/tensor.h"

namespace fluoroseq {

DetachTransition::DetachTransition(double p_detach) : p_detach(p_detach) {}

void DetachTransition::operator()(Tensor* tensor, int timestep) const {
    for (int t = 0; t <= timestep; t++) {
        double sum = 0.0;
        for (int i = 1; i < tensor->strides[0]; i++) {
            double value = tensor->values[t * tensor->strides[0] + i];
            tensor->values[t * tensor->strides[0] + i] = value * (1 - p_detach);
            sum += value;
        }
        tensor->values[t * tensor->strides[0]] += p_detach * sum;
    }
}

}  // namespace fluoroseq
