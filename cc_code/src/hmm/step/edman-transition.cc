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
#include "hmm/state-vector/peptide-state-vector.h"
#include "parameterization/fit/sequencing-model-fitter.h"
#include "tensor/vector.h"
#include "util/kd-range.h"

namespace whatprot {

EdmanTransition::EdmanTransition(double p_edman_failure,
                                 const DyeSeq& dye_seq,
                                 const DyeTrack& dye_track)
        : dye_seq(dye_seq),
          dye_track(dye_track),
          p_edman_failure(p_edman_failure) {}

void EdmanTransition::prune_forward(KDRange* range, bool* allow_detached) {
    forward_range = *range;
    range->max[0]++;
    for (unsigned int c = 0; c < range->min.size() - 1; c++) {
        if (range->min[1 + c] != 0) {
            range->min[1 + c]--;
        }
    }
    backward_range = *range;
}

void EdmanTransition::prune_backward(KDRange* range, bool* allow_detached) {
    backward_range = backward_range.intersect(*range);
    *range = backward_range;
    if (range->min[0] != 0) {
        range->min[0]--;
    }
    for (unsigned int c = 0; c < range->min.size() - 1; c++) {
        range->max[1 + c]++;
    }
    *range = forward_range.intersect(*range);
    forward_range = *range;
}

PeptideStateVector* EdmanTransition::forward(const PeptideStateVector& input,
                                             unsigned int* num_edmans) const {
    (*num_edmans)++;
    PeptideStateVector* output = new PeptideStateVector(
            backward_range.max.size(), &backward_range.max[0]);
    // First we set all of the output in the backward range to zero. This allows
    // us to use += when gathering the various probabilities coming in from the
    // input PeptideStateVector. This is way easier than the alternative, since
    // the forward and backward ranges may not match up the way you expect.
    TensorIterator* out_itr = output->tensor.iterator(backward_range);
    while (!out_itr->done()) {
        *out_itr->get() = 0.0;
        out_itr->advance();
    }
    delete out_itr;
    // Now we iterate through the input in the forward range, and multiply these
    // values out into every receiving value in the output. We don't worry about
    // whether we are writing to locations which are actually in the
    // backward range, as this would be more trouble than it's worth, and likely
    // would not improve the runtime (checking conditionals is expensive).
    ConstTensorIterator* in_itr = input.tensor.const_iterator(forward_range);
    unsigned int t_stride = output->tensor.strides[0];
    while (!in_itr->done()) {
        double f_val = *in_itr->get();
        unsigned int i = in_itr->index;
        unsigned int t = in_itr->loc[0];
        int c = dye_seq[t];
        // Probability of failure is straightforward.
        output->tensor.values[i] += p_edman_failure * f_val;
        // Probability of success is broken into smaller pieces.
        if (c == -1) {
            // If no fluorophore removed in the successful Edman cycle scenario,
            // we just take the remaining portion of the probability to the next
            // successful Edman count ('+ t_stride' does this)
            output->tensor.values[i + t_stride] +=
                    (1 - p_edman_failure) * f_val;
        } else {
            // If fluorophore removed, we need to split this probability
            // further. As before we reindex to the next successful Edman count
            // (with '+ t_stride').
            unsigned int c_idx = in_itr->loc[1 + c];
            unsigned int c_total = dye_track(t, c);
            double ratio = (double)c_idx / (double)c_total;
            unsigned int c_stride = output->tensor.strides[1 + c];
            if (c_idx < c_total) {
                // Here we multiply additionally by probability of no
                // fluorophore removal.
                output->tensor.values[i + t_stride] +=
                        (1 - p_edman_failure) * (1 - ratio) * f_val;
            }
            if (c_idx > 0) {
                // And here we handle probability of flurophore removal
                // ('- c_stride' indexes to the correct location).
                output->tensor.values[i + t_stride - c_stride] +=
                        (1 - p_edman_failure) * ratio * f_val;
            }
        }
        in_itr->advance();
    }
    delete in_itr;
    output->range = backward_range;
    output->allow_detached = input.allow_detached;
    if (output->allow_detached) {
        output->p_detached = input.p_detached;
    }
    return output;
}

PeptideStateVector* EdmanTransition::backward(const PeptideStateVector& input,
                                              unsigned int* num_edmans) const {
    PeptideStateVector* output = new PeptideStateVector(
            backward_range.max.size(), &backward_range.max[0]);
    // First we set all of the output in the forward range to zero. This allows
    // us to use += when gathering the various probabilities coming in from the
    // input PeptideStateVector. This is way easier than the alternative, since
    // the forward and backward ranges may not match up the way you expect.
    TensorIterator* out_itr = output->tensor.iterator(forward_range);
    while (!out_itr->done()) {
        *out_itr->get() = 0.0;
        out_itr->advance();
    }
    delete out_itr;
    // Now we iterate through the input in the backward range, and multiply
    // these values out into every receiving value in the output. We don't worry
    // about whether we are writing to locations which are actually in the
    // backward range, as this would be more trouble than it's worth, and likely
    // would not improve the runtime (checking conditionals is expensive).
    ConstTensorIterator* in_itr = input.tensor.const_iterator(backward_range);
    unsigned int t_stride = output->tensor.strides[0];
    while (!in_itr->done()) {
        double f_val = *in_itr->get();
        unsigned int i = in_itr->index;
        unsigned int t = in_itr->loc[0];
        // Probability of failure is straightforward.
        output->tensor.values[i] += p_edman_failure * f_val;
        // Probability of success is broken into smaller pieces. Note that this
        // is only relevant when t > 0, because the previous number of
        // successful Edman cycles is t - 1, which is invalid at -1.
        if (t > 0) {
            int c = dye_seq[t - 1];
            if (c == -1) {
                // If no fluorophore removed in the successful Edman cycle
                // scenario, we just take the remaining portion of the
                // probability to the next successful Edman count ('- t_stride'
                // does this).
                output->tensor.values[i - t_stride] +=
                        (1 - p_edman_failure) * f_val;
            } else {
                // If fluorophore removed, we need to split this probability
                // further. As before we reindex to the next successful Edman
                // count (with '- t_stride'). Note that we get different values
                // for 'c_idx' and 'ratio' in the two cases, because we are
                // computing 'backward' based on 'forward'. The 'forward' part
                // may either have the same number of fluorophores as the target
                // location in 'backward' (thus not losing a dye) or one more
                // fluorophore, which it loses.
                unsigned int c_total = dye_track(t - 1, c);
                unsigned int c_stride = output->tensor.strides[1 + c];
                unsigned int c_idx;
                c_idx = in_itr->loc[1 + c];
                if (c_idx < c_total) {
                    // Here we multiply additionally by probability of no
                    // fluorophore removal.
                    double ratio = (double)c_idx / (double)c_total;
                    output->tensor.values[i - t_stride] +=
                            (1 - p_edman_failure) * (1 - ratio) * f_val;
                }
                c_idx = in_itr->loc[1 + c] + 1;
                if (c_idx > 0) {
                    // And here we handle probability of flurophore removal
                    // ('+ c_stride' indexes to the correct location).
                    double ratio = (double)c_idx / (double)c_total;
                    output->tensor.values[i - t_stride + c_stride] +=
                            (1 - p_edman_failure) * ratio * f_val;
                }
            }
        }
        in_itr->advance();
    }
    delete in_itr;
    output->range = forward_range;
    output->allow_detached = input.allow_detached;
    if (output->allow_detached) {
        output->p_detached = input.p_detached;
    }
    (*num_edmans)--;
    return output;
}

void EdmanTransition::improve_fit(const PeptideStateVector& forward_psv,
                                  const PeptideStateVector& backward_psv,
                                  const PeptideStateVector& next_backward_psv,
                                  unsigned int num_edmans,
                                  double probability,
                                  SequencingModelFitter* fitter) const {
    unsigned int t_stride = forward_psv.tensor.strides[0];
    for (unsigned int t = 0; t < num_edmans + 1; t++) {
        // Here we omit the zeroth entry of every timestep because this is the
        // entry for zero of every dye color. These entries are unable to
        // provide tangible evidence of the edman efficiency one way or the
        // other.
        for (unsigned int i = t * t_stride + 1; i < (t + 1) * t_stride; i++) {
            fitter->p_edman_failure_fit.numerator +=
                    forward_psv.tensor.values[i] * p_edman_failure
                    * next_backward_psv.tensor.values[i] / probability;
            fitter->p_edman_failure_fit.denominator +=
                    forward_psv.tensor.values[i] * backward_psv.tensor.values[i]
                    / probability;
        }
    }
}

}  // namespace whatprot
