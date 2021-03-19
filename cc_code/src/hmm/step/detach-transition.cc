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
#include "hmm/state-vector/peptide-state-vector.h"

namespace whatprot {

DetachTransition::DetachTransition(double p_detach) : p_detach(p_detach) {}

void DetachTransition::forward(PeptideStateVector* psv) const {
    int i_max = (psv->num_edmans + 1) * psv->tensor.strides[0];
    double sum = 0.0;
    for (int i = 0; i < i_max; i++) {
        double value = psv->tensor.values[i];
        psv->tensor.values[i] = value * (1 - p_detach);
        sum += value;
    }
    psv->tensor.values[psv->num_edmans * psv->tensor.strides[0]] +=
            p_detach * sum;
}

void DetachTransition::backward(const PeptideStateVector& input,
                                PeptideStateVector* output) const {
    int i_max = (input.num_edmans + 1) * input.tensor.strides[0];
    double if_detach =
            input.tensor.values[input.num_edmans * input.tensor.strides[0]];
    for (int i = 0; i < i_max; i++) {
        output->tensor.values[i] =
                (1 - p_detach) * input.tensor.values[i] + p_detach * if_detach;
    }
    output->num_edmans = input.num_edmans;
}

void DetachTransition::improve_fit(const PeptideStateVector& forward_psv,
                                   const PeptideStateVector& backward_psv,
                                   const PeptideStateVector& next_backward_psv,
                                   double probability,
                                   ErrorModelFitter* fitter) const {
    int t_stride = forward_psv.tensor.strides[0];
    double forward_sum = 0.0;
    double forward_backward_sum = 0.0;
    for (int t = 0; t < forward_psv.num_edmans + 1; t++) {
        // Here we omit the zeroth entry of every timestep because this is the
        // entry for zero of every dye color. These entries are unable to
        // provide tangible evidence of detachment one way or the other.
        for (int i = t * t_stride + 1; i < (t + 1) * t_stride; i++) {
            forward_sum += forward_psv.tensor.values[i];
            forward_backward_sum += forward_psv.tensor.values[i]
                                    * backward_psv.tensor.values[i];
        }
    }
    fitter->p_detach_fit.numerator +=
            forward_sum * p_detach
            * next_backward_psv.tensor.values[forward_psv.num_edmans
                                              * forward_psv.tensor.strides[0]]
            / probability;
    // Probability of being in a state that can detach is 1.0, because all
    // states can detach (although we are ignoring the case where there are no
    // amino acids left but this probably shouldn't cause any serious issues).
    fitter->p_detach_fit.denominator += forward_backward_sum / probability;
}

}  // namespace whatprot
