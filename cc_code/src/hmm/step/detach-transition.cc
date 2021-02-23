/******************************************************************************\
* Author: Matthew Beauregard Smith                                             *
* Affiliation: The University of Texas at Austin                               *
* Department: Oden Institute and Institute for Cellular and Molecular Biology  *
* PI: Edward Marcotte                                                          *
* Project: Protein Fluorosequencing                                            *
\******************************************************************************/

// Defining symbols from header:
#include "detach-transition.h"

// Local project headers:
#include "hmm/fit/error-model-fitter.h"
#include "tensor/tensor.h"

namespace whatprot {

DetachTransition::DetachTransition(double p_detach) : p_detach(p_detach) {}

void DetachTransition::forward(int* edmans,
                               Tensor* output) const {
    int i_max = ((*edmans) + 1) * output->strides[0];
    double sum = 0.0;
    for (int i = 0; i < i_max; i++) {
        double value = output->values[i];
        output->values[i] = value * (1 - p_detach);
        sum += value;
    }
    output->values[(*edmans) * output->strides[0]] += p_detach * sum;
}

void DetachTransition::backward(const Tensor& input,
                                int* edmans,
                                Tensor* output) const {
    int i_max = ((*edmans) + 1) * input.strides[0];
    double if_detach = input.values[(*edmans) * input.strides[0]];
    for (int i = 0; i < i_max; i++) {
        output->values[i] =
                (1 - p_detach) * input.values[i] + p_detach * if_detach;
    }
}

void DetachTransition::improve_fit(const Tensor& forward_tensor,
                                   const Tensor& backward_tensor,
                                   const Tensor& next_backward_tensor,
                                   int edmans,
                                   double probability,
                                   ErrorModelFitter* fitter) const {
    int t_stride = forward_tensor.strides[0];
    double forward_sum = 0.0;
    double forward_backward_sum = 0.0;
    for (int t = 0; t < edmans + 1; t++) {
        // Here we omit the zeroth entry of every timestep because this is the
        // entry for zero of every dye color. These entries are unable to
        // provide tangible evidence of detachment one way or the other.
        for (int i = t * t_stride + 1; i < (t + 1) * t_stride; i++) {
            forward_sum += forward_tensor.values[i];
            forward_backward_sum +=
                    forward_tensor.values[i] * backward_tensor.values[i];
        }
    }
    fitter->p_detach_fit.numerator +=
            forward_sum * p_detach
            * next_backward_tensor.values[edmans * forward_tensor.strides[0]]
            / probability;
    // Probability of being in a state that can detach is 1.0, because all
    // states can detach (although we are ignoring the case where there are no
    // amino acids left but this probably shouldn't cause any serious issues).
    fitter->p_detach_fit.denominator += forward_backward_sum / probability;
}

}  // namespace whatprot
