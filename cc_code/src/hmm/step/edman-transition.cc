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

void EdmanTransition::set_true_forward_range(const KDRange& range) {
    true_forward_range = range;
    safe_backward_range = range;
    safe_backward_range.max[0]++;
    for (unsigned int c = 0; c < safe_backward_range.min.size() - 1; c++) {
        if (safe_backward_range.min[1 + c] != 0) {
            safe_backward_range.min[1 + c]--;
        }
    }
}

void EdmanTransition::set_true_backward_range(const KDRange& range) {
    true_backward_range = range;
    safe_forward_range = range;
    if (safe_forward_range.min[0] != 0) {
        safe_forward_range.min[0]--;
    }
    for (unsigned int c = 0; c < safe_forward_range.min.size() - 1; c++) {
        safe_forward_range.max[1 + c]++;
    }
}

void EdmanTransition::prune_forward(KDRange* range, bool* allow_detached) {
    set_true_forward_range(*range);
    *range = safe_backward_range;
}

void EdmanTransition::prune_backward(KDRange* range, bool* allow_detached) {
    *range = safe_backward_range.intersect(*range);
    set_true_backward_range(*range);
    *range = safe_forward_range.intersect(true_forward_range);
    set_true_forward_range(*range);
}

PeptideStateVector* EdmanTransition::forward(const PeptideStateVector& input,
                                             unsigned int* num_edmans) const {
    (*num_edmans)++;
    PeptideStateVector* output = new PeptideStateVector(safe_backward_range);
    // First we set all of the output in the backward range to zero. This allows
    // us to use += when gathering the various probabilities coming in from the
    // input PeptideStateVector. This is way easier than the alternative, since
    // the forward and backward ranges may not match up the way you expect.
    TensorIterator* out_itr = output->tensor.iterator(safe_backward_range);
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
    ConstTensorIterator* in_itr =
            input.tensor.const_iterator(true_forward_range);
    // true_forward_range is a strict subset of safe_backward_range, so we can
    // use it to index into output.
    out_itr = output->tensor.iterator(true_forward_range);
    unsigned int t_stride = output->tensor.strides[0];
    while (!in_itr->done()) {
        double f_val = *in_itr->get();
        unsigned int i = out_itr->index;
        unsigned int t = out_itr->loc[0];
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
        out_itr->advance();
    }
    delete in_itr;
    delete out_itr;
    // We also need to deal with the 'broken-n' states. Even though Edman
    // degradation has no effect on these states, which is the whole reason for
    // their existence, we still need to copy them to the new tensor so that the
    // old values are not lost.
    //
    // Note that just as with the normal tensor states, we first just set every
    // value in the 'safe-backward-range' to zero, as that is much easier than
    // tracking everything properly.
    TensorIterator* out_n_itr =
            output->broken_n_tensor.iterator(safe_backward_range);
    while (!out_n_itr->done()) {
        *out_n_itr->get() = 0.0;
        out_n_itr->advance();
    }
    delete out_n_itr;
    // Now we actually transfer the values.
    ConstTensorIterator* in_n_itr =
            input.broken_n_tensor.const_iterator(true_forward_range);
    // true_forward_range is a strict subset of safe_backward_range, so we can
    // use it to index into output.
    out_n_itr = output->broken_n_tensor.iterator(true_forward_range);
    while (!in_n_itr->done()) {
        *out_n_itr->get() = *in_n_itr->get();
        in_n_itr->advance();
        out_n_itr->advance();
    }
    delete in_n_itr;
    delete out_n_itr;
    // Now we fix up the ranges, allow_detached, etc...
    output->range = true_backward_range;
    output->allow_detached = input.allow_detached;
    if (output->allow_detached) {
        output->p_detached = input.p_detached;
    }
    return output;
}

PeptideStateVector* EdmanTransition::backward(const PeptideStateVector& input,
                                              unsigned int* num_edmans) const {
    PeptideStateVector* output = new PeptideStateVector(safe_forward_range);
    // First we set all of the output in the forward range to zero. This allows
    // us to use += when gathering the various probabilities coming in from the
    // input PeptideStateVector. This is way easier than the alternative, since
    // the forward and backward ranges may not match up the way you expect.
    TensorIterator* out_itr = output->tensor.iterator(safe_forward_range);
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
    ConstTensorIterator* in_itr =
            input.tensor.const_iterator(true_backward_range);
    // true_backward_range is a strict subset of safe_forward_range, so we can
    // use it to index into output.
    out_itr = output->tensor.iterator(true_backward_range);
    unsigned int t_stride = output->tensor.strides[0];
    while (!in_itr->done()) {
        double f_val = *in_itr->get();
        unsigned int i = out_itr->index;
        unsigned int t = out_itr->loc[0];
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
        out_itr->advance();
    }
    delete in_itr;
    delete out_itr;
    // We also need to deal with the 'broken-n' states. Even though Edman
    // degradation has no effect on these states, which is the whole reason for
    // their existence, we still need to copy them to the new tensor so that the
    // old values are not lost.
    //
    // Note that just as with the normal tensor states, we first just set every
    // value in the 'safe-backward-range' to zero, as that is much easier than
    // tracking everything properly.
    TensorIterator* out_n_itr =
            output->broken_n_tensor.iterator(safe_forward_range);
    while (!out_n_itr->done()) {
        *out_n_itr->get() = 0.0;
        out_n_itr->advance();
    }
    delete out_n_itr;
    // Now we actually transfer the values.
    ConstTensorIterator* in_n_itr =
            input.broken_n_tensor.const_iterator(true_backward_range);
    // true_backward_range is a strict subset of safe_forward_range, so we can
    // use it to index into output.
    out_n_itr = output->broken_n_tensor.iterator(true_backward_range);
    while (!in_n_itr->done()) {
        *out_n_itr->get() = *in_n_itr->get();
        in_n_itr->advance();
        out_n_itr->advance();
    }
    delete in_n_itr;
    delete out_n_itr;
    // Now we fix up the ranges, allow_detached, etc...
    output->range = true_forward_range;
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
    ConstTensorIterator* f_itr =
            forward_psv.tensor.const_iterator(true_forward_range);
    ConstTensorIterator* b_itr =
            backward_psv.tensor.const_iterator(true_forward_range);
    ConstTensorIterator* nb_itr =
            next_backward_psv.tensor.const_iterator(true_forward_range);
    unsigned int old_t = -1;
    while (!f_itr->done()) {
        // Here we omit the zeroth entry of every timestep because this is the
        // entry for zero of every dye color. These entries are unable to
        // provide tangible evidence of the edman efficiency one way or the
        // other.
        if (true_forward_range.includes_zero()) {
            unsigned int new_t = f_itr->loc[0];
            if (new_t != old_t) {
                old_t = new_t;
                f_itr->advance();
                b_itr->advance();
                nb_itr->advance();
                // Now we need to start the loop from the beginning; maybe all
                // entries are zero dyes on all channels, or maybe we have met
                // the loop end condition. Either way we need to check.
                continue;
            }
        }
        // And now we can accumulate information about Edman failure rate.
        double f_prob = *f_itr->get();
        double b_prob = *b_itr->get();
        double nb_prob = *nb_itr->get();
        fitter->p_edman_failure_fit.numerator +=
                f_prob * p_edman_failure * nb_prob / probability;
        fitter->p_edman_failure_fit.denominator +=
                f_prob * b_prob / probability;
        f_itr->advance();
        b_itr->advance();
        nb_itr->advance();
    }
    delete f_itr;
    delete b_itr;
    delete nb_itr;
}

}  // namespace whatprot
