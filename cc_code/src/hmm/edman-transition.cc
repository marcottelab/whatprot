/******************************************************************************\
* Author: Matthew Beauregard Smith                                             *
* Affiliation: The University of Texas at Austin                               *
* Department: Oden Institute and Institute for Cellular and Molecular Biology  *
* PI: Edward Marcotte                                                          *
* Project: Protein Fluorosequencing                                            *
\******************************************************************************/

// Defining symbols from header:
#include "edman-transition.h"

// Local project headers:
#include "common/dye-track.h"
#include "tensor/tensor.h"
#include "tensor/vector.h"

namespace fluoroseq {

EdmanTransition::EdmanTransition(double p_edman_failure,
                                 const DyeSeq& dye_seq,
                                 const DyeTrack& dye_track)
        : p_edman_failure(p_edman_failure),
          dye_seq(dye_seq),
          dye_track(dye_track) {}

void EdmanTransition::forward(const Tensor& input, int timestep, Tensor* output) const {
    int t_stride = input.strides[0];
    for (int i = t_stride * (1 + timestep) - 1; i >= 0; i--) {
        output->values[i + t_stride] = input.values[i];
    }
    for (int i = 0; i < t_stride; i++) {
        output->values[i] *= p_edman_failure;
    }
    for (int t = 0; t < timestep + 1; t++) {
        if (t > 0) {
            for (int i = t * t_stride; i < (t + 1) * t_stride; i++) {
                output->values[i] +=
                        p_edman_failure * input.values[i + t_stride];
            }
        }
        short channel = dye_seq[t];
        if (channel != -1) {
            int amt = dye_track(t, channel);
            int vector_stride = input.strides[1 + channel];
            int vector_length = input.shape[1 + channel];
            int outer_stride = vector_stride * vector_length;
            int outer_min = (t + 1) * t_stride;
            int outer_max = (t + 2) * t_stride;
            for (int outer = outer_min; outer < outer_max;
                 outer += outer_stride) {
                for (int inner = 0; inner < vector_stride; inner++) {
                    const Vector inv(vector_length,
                             vector_stride,
                             &input.values[outer + inner]);
                    Vector outv(vector_length,
                            vector_stride,
                            &output->values[outer + inner]);
                    for (int i = 1; i <= amt; i++) {
                        double ratio = (double)i / (double)amt;
                        outv[i - 1] += inv[i] * ratio;
                        outv[i] *= 1 - ratio;
                    }
                }
            }
        }
        for (int i = (t + 1) * t_stride; i < (t + 2) * t_stride; i++) {
            output->values[i] *= (1 - p_edman_failure);
        }
    }
}

}  // namespace fluoroseq
