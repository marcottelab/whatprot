// Author: Matthew Beauregard Smith
#include "edman_transition.h"

#include "common/dye_track.h"
#include "tensor/tensor.h"
#include "tensor/vector.h"

namespace fluoroseq {

EdmanTransition::EdmanTransition(double p_edman_failure,
                                 const DyeSeq& dye_seq,
                                 const DyeTrack& dye_track)
        : p_edman_failure(p_edman_failure),
          dye_seq(dye_seq),
          dye_track(dye_track) {}

void EdmanTransition::operator()(Tensor* tensor, int timestep) const {
    int t_stride = tensor->strides[0];
    for (int i = t_stride * (1 + timestep) - 1; i >= 0; i--) {
        tensor->values[i + t_stride] = tensor->values[i];
    }
    for (int i = 0; i < t_stride; i++) {
        tensor->values[i] *= p_edman_failure;
    }
    for (int t = 0; t < timestep + 1; t++) {
        if (t > 0) {
            for (int i = t * t_stride; i < (t + 1) * t_stride; i++) {
                tensor->values[i] += p_edman_failure
                                     * tensor->values[i + t_stride];
            }
        }
        short channel = dye_seq[t];
        if (channel != -1) {
            int amt = dye_track(t, channel);
            int vector_stride = tensor->strides[1 + channel];
            int vector_length = tensor->shape[1 + channel];
            int outer_stride = vector_stride * vector_length;
            int outer_min = (t + 1) * t_stride;
            int outer_max = (t + 2) * t_stride;
            for (int outer = outer_min;
                     outer < outer_max;
                     outer += outer_stride) {
                for (int inner = 0; inner < vector_stride; inner++) {
                    Vector v(vector_length,
                             vector_stride,
                             &tensor->values[outer+inner]);
                    for (int i = 1; i <= amt; i++) {
                        double ratio = (double) i / (double) amt;
                        v[i - 1] += v[i] * ratio;
                        v[i] *= 1 - ratio;
                    }
                }
            }
        }
        for (int i = (t + 1) * t_stride; i < (t + 2) * t_stride; i++) {
            tensor->values[i] *= (1 - p_edman_failure);
        }
    }
}

}  // namespace fluoroseq
