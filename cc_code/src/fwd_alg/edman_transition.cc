// Author: Matthew Beauregard Smith
#include "edman_transition.h"

#include <algorithm>

#include "common/dye_track.h"
#include "tensor/tensor.h"
#include "tensor/vector.h"

namespace fluoroseq {

namespace {
using std::copy;
}  // namespace

EdmanTransition::EdmanTransition(double p_edman_failure,
                                 const DyeSeq& dye_seq,
                                 const DyeTrack& dye_track,
                                 int max_failed_edmans)
        : p_edman_failure(p_edman_failure),
          dye_seq(dye_seq),
          dye_track(dye_track),
          max_failed_edmans(max_failed_edmans) {}

void EdmanTransition::operator()(Tensor* tensor, int timestep) const {
    int t_stride = tensor->strides[0];
    int t_min = timestep - max_failed_edmans - 1;
    bool removing_lowest = (t_min >= 0);
    if (t_min < 0) {
        t_min = 0;
    }
    for (int i = t_stride * timestep - 1; i >= t_min * t_stride; i--) {
        tensor->values[i + t_stride] = tensor->values[i];
    }
    if (!removing_lowest) {
        for (int i = 0; i < t_stride; i++) {
            tensor->values[i] *= p_edman_failure;
        }
    }
    for (int t = t_min; t < timestep; t++) {
        if (t > t_min) {
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
