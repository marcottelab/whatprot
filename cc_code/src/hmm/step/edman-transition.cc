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
#include "hmm/fit/error-model-fitter.h"
#include "tensor/tensor.h"
#include "tensor/vector.h"

namespace whatprot {

EdmanTransition::EdmanTransition(double p_edman_failure,
                                 const DyeSeq& dye_seq,
                                 const DyeTrack& dye_track)
        : p_edman_failure(p_edman_failure),
          dye_seq(dye_seq),
          dye_track(dye_track) {}

void EdmanTransition::forward(int* edmans,
                              Tensor* output) const {
    (*edmans)++;
    int t_stride = output->strides[0];
    for (int i = t_stride * (*edmans) - 1; i >= 0; i--) {
        output->values[i + t_stride] = output->values[i];
    }
    for (int i = 0; i < t_stride; i++) {
        output->values[i] = output->values[i] * p_edman_failure;
    }
    for (int t = 0; t < *edmans; t++) {
        if (t > 0) {
            for (int i = t * t_stride; i < (t + 1) * t_stride; i++) {
                output->values[i] +=
                        p_edman_failure * output->values[i + t_stride];
            }
        }
        short channel = dye_seq[t];
        if (channel != -1) {
            int amt = dye_track(t, channel);
            int vector_stride = output->strides[1 + channel];
            int vector_length = output->shape[1 + channel];
            int outer_stride = vector_stride * vector_length;
            int outer_min = (t + 1) * t_stride;
            int outer_max = (t + 2) * t_stride;
            for (int outer = outer_min; outer < outer_max;
                 outer += outer_stride) {
                for (int inner = 0; inner < vector_stride; inner++) {
                    Vector v(vector_length,
                             vector_stride,
                             &output->values[outer + inner]);
                    for (int i = 1; i <= amt; i++) {
                        double ratio = (double)i / (double)amt;
                        v[i - 1] += v[i] * ratio;
                        v[i] *= 1 - ratio;
                    }
                }
            }
        }
        for (int i = (t + 1) * t_stride; i < (t + 2) * t_stride; i++) {
            output->values[i] *= (1 - p_edman_failure);
        }
    }
}

void EdmanTransition::backward(const Tensor& input,
                               int* edmans,
                               Tensor* output) const {
    int t_stride = input.strides[0];
    for (int t = 0; t < *edmans; t++) {
        for (int i = t * t_stride; i < (t + 1) * t_stride; i++) {
            output->values[i] = p_edman_failure * input.values[i];
        }
        short channel = dye_seq[t];
        if (channel == -1) {
            for (int i = t * t_stride; i < (t + 1) * t_stride; i++) {
                output->values[i] +=
                        (1 - p_edman_failure) * input.values[i + t_stride];
            }
        } else {
            int amt = dye_track(t, channel);
            int vector_stride = output->strides[1 + channel];
            int vector_length = output->shape[1 + channel];
            int outer_stride = vector_stride * vector_length;
            int outer_min = t * t_stride;
            int outer_max = (t + 1) * t_stride;
            for (int outer = outer_min; outer < outer_max;
                 outer += outer_stride) {
                for (int inner = 0; inner < vector_stride; inner++) {
                    const Vector inv(vector_length,
                                     vector_stride,
                                     &input.values[outer + inner + t_stride]);
                    Vector outv(vector_length,
                                vector_stride,
                                &output->values[outer + inner]);
                    for (int i = 0; i <= amt; i++) {
                        double ratio = (double)i / (double)amt;
                        if (i != 0) {
                            outv[i] +=
                                    (1 - p_edman_failure) * ratio * inv[i - 1];
                        }
                        outv[i] += (1 - p_edman_failure) * (1 - ratio) * inv[i];
                    }
                }
            }
        }
    }
    (*edmans)--;
}

void EdmanTransition::improve_fit(const Tensor& forward_tensor,
                                  const Tensor& backward_tensor,
                                  const Tensor& next_backward_tensor,
                                  int edmans,
                                  double probability,
                                  ErrorModelFitter* fitter) const {
    int t_stride = forward_tensor.strides[0];
    for (int t = 0; t < edmans + 1; t++) {
        // Here we omit the zeroth entry of every timestep because this is the
        // entry for zero of every dye color. These entries are unable to
        // provide tangible evidence of the edman efficiency one way or the
        // other.
        for (int i = t * t_stride + 1; i < (t + 1) * t_stride; i++) {
            fitter->p_edman_failure_fit.numerator +=
                    forward_tensor.values[i] * p_edman_failure
                    * next_backward_tensor.values[i] / probability;
            fitter->p_edman_failure_fit.denominator +=
                    forward_tensor.values[i] * backward_tensor.values[i]
                    / probability;
        }
    }
}

}  // namespace whatprot
