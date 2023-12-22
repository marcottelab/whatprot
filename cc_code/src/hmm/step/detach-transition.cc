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
#include "hmm/state-vector/peptide-state-vector.h"
#include "parameterization/fit/sequencing-model-fitter.h"
#include "tensor/const-tensor-iterator.h"
#include "tensor/tensor-iterator.h"
#include "util/kd-range.h"

namespace whatprot {

DetachTransition::DetachTransition(unsigned int timestep, double p_detach)
        : p_detach(p_detach), timestep(timestep) {}

void DetachTransition::prune_forward(KDRange* range, bool* allow_detached) {
    pruned_range = *range;
    detached_forward = *allow_detached;
    *allow_detached = true;
}

void DetachTransition::prune_backward(KDRange* range, bool* allow_detached) {
    pruned_range = pruned_range.intersect(*range);
    *range = pruned_range;
    detached_backward = allow_detached;
}

PeptideStateVector* DetachTransition::forward(const PeptideStateVector& input,
                                              unsigned int* num_edmans) const {
    PeptideStateVector* output = new PeptideStateVector(pruned_range);
    double sum = forward(input.tensor, &output->tensor);
    sum += forward(input.broken_n_tensor, &output->broken_n_tensor);
    if (detached_backward) {
        if (detached_forward) {
            output->p_detached = input.p_detached + p_detach * sum;
        } else {
            output->p_detached = p_detach * sum;
        }
    }
    // Now we fix up the ranges, allow_detached, etc...
    output->range = pruned_range;
    output->allow_detached = detached_backward;
    return output;
}

double DetachTransition::forward(const Tensor& input, Tensor* output) const {
    ConstTensorIterator* in_itr = input.const_iterator(pruned_range);
    TensorIterator* out_itr = output->iterator(pruned_range);
    double sum = 0.0;
    while (!in_itr->done()) {
        double value = *in_itr->get();
        *out_itr->get() = value * (1 - p_detach);
        sum += value;
        in_itr->advance();
        out_itr->advance();
    }
    delete in_itr;
    delete out_itr;
    return sum;
}

PeptideStateVector* DetachTransition::backward(const PeptideStateVector& input,
                                               unsigned int* num_edmans) const {
    PeptideStateVector* output = new PeptideStateVector(pruned_range);
    backward(input.tensor, input.p_detached, &output->tensor);
    backward(input.broken_n_tensor, input.p_detached, &output->broken_n_tensor);
    if (detached_forward) {
        if (detached_backward) {
            output->p_detached = input.p_detached;
        } else {
            output->p_detached = 0.0;
        }
    }
    // Now we fix up the ranges, allow_detached, etc...
    output->range = pruned_range;
    output->allow_detached = detached_forward;
    return output;
}

void DetachTransition::backward(const Tensor& input,
                                double p_detached,
                                Tensor* output) const {
    ConstTensorIterator* in_itr = input.const_iterator(pruned_range);
    TensorIterator* out_itr = output->iterator(pruned_range);
    while (!in_itr->done()) {
        *out_itr->get() = (1 - p_detach) * (*in_itr->get());
        if (detached_backward) {
            *out_itr->get() += p_detach * p_detached;
        }
        in_itr->advance();
        out_itr->advance();
    }
    delete in_itr;
    delete out_itr;
}

void DetachTransition::improve_fit(const PeptideStateVector& forward_psv,
                                   const PeptideStateVector& backward_psv,
                                   const PeptideStateVector& next_backward_psv,
                                   unsigned int num_edmans,
                                   double probability,
                                   SequencingModelFitter* fitter) const {
    improve_fit(forward_psv.tensor,
                backward_psv.tensor,
                next_backward_psv.p_detached,
                probability,
                fitter);
    improve_fit(forward_psv.broken_n_tensor,
                backward_psv.broken_n_tensor,
                next_backward_psv.p_detached,
                probability,
                fitter);
}

void DetachTransition::improve_fit(const Tensor& forward_tsr,
                                   const Tensor& backward_tsr,
                                   double next_backward_p_detached,
                                   double probability,
                                   SequencingModelFitter* fitter) const {
    ConstTensorIterator* f_itr = forward_tsr.const_iterator(pruned_range);
    ConstTensorIterator* b_itr = backward_tsr.const_iterator(pruned_range);
    double forward_sum = 0.0;
    double forward_backward_sum = 0.0;
    unsigned int old_t = -1;
    while (!f_itr->done()) {
        // Here we omit the zeroth entry of every timestep because this is the
        // entry for zero of every dye color. These entries are unable to
        // provide tangible evidence of detachment one way or the other.
        if (pruned_range.includes_zero()) {
            unsigned int new_t = f_itr->loc[0];
            if (new_t != old_t) {
                old_t = new_t;
                f_itr->advance();
                b_itr->advance();
                // Now we need to start the loop from the beginning; maybe all
                // entries are zero dyes on all channels, or maybe we have met
                // the loop end condition. Either way we need to check.
                continue;
            }
        }
        // And now we can accumulate forward and forward-backward sums
        double f_prob = *f_itr->get();
        double b_prob = *b_itr->get();
        forward_sum += f_prob;
        forward_backward_sum += f_prob * b_prob;
        f_itr->advance();
        b_itr->advance();
    }
    delete f_itr;
    delete b_itr;
    double numerator =
            forward_sum * p_detach * next_backward_p_detached / probability;
    // Probability of being in a state that can detach is 1.0, because all
    // states can detach (although we are ignoring the case where there are no
    // amino acids left but that's fine and in fact better).
    double denominator = forward_backward_sum / probability;
    fitter->p_detach_fit.add_timestep(timestep, numerator, denominator);
}

}  // namespace whatprot
