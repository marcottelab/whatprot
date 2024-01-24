/******************************************************************************\
* Author: Matthew Beauregard Smith                                             *
* Affiliation: The University of Texas at Austin                               *
* Department: Oden Institute and Institute for Cellular and Molecular Biology  *
* PI: Edward Marcotte                                                          *
* Project: Protein Fluorosequencing                                            *
\******************************************************************************/

// Defining symbols from header:
#include "peptide-emission.h"

// Standard C++ library headers:
#include <algorithm>
#include <cmath>
#include <limits>
#include <vector>

// Local project headers:
#include "common/radiometry.h"
#include "hmm/state-vector/peptide-state-vector.h"
#include "parameterization/fit/sequencing-model-fitter.h"
#include "parameterization/model/channel-model.h"
#include "parameterization/model/sequencing-model.h"
#include "parameterization/settings/sequencing-settings.h"
#include "tensor/const-tensor-iterator.h"
#include "tensor/tensor-iterator.h"
#include "util/kd-range.h"

namespace whatprot {

namespace {
using std::function;
using std::lround;
using std::min;
using std::numeric_limits;
using std::vector;
}  // namespace

PeptideEmission::PeptideEmission(const Radiometry& radiometry,
                                 unsigned int timestep,
                                 unsigned int max_num_dyes,
                                 const SequencingModel& seq_model,
                                 const SequencingSettings& seq_settings)
        : radiometry(radiometry),
          timestep(timestep),
          i_am_a_copy(false),
          num_channels(radiometry.num_channels),
          max_num_dyes(max_num_dyes) {
    pruned_range.min.resize(1 + num_channels);
    pruned_range.max.resize(1 + num_channels);
    pruned_range.min[0] = 0;
    pruned_range.max[0] = timestep + 1;
    if (seq_settings.dist_cutoff == numeric_limits<double>::max()) {
        for (unsigned int c = 0; c < num_channels; c++) {
            pruned_range.min[1 + c] = 0;
            pruned_range.max[1 + c] = max_num_dyes + 1;
        }
    } else {
        // This vector is a rough estimate of what the lower-bound of the
        // pruned_range will be for channel c. This is used to call the
        // adjusted_mu() function to determine mu for a channel, given FRET
        // effects with other channels and quenching effects with itself. We
        // just set this one to all zeros. Something similar was tried to what
        // is done for counts_hi (see below), but this seemed to result in an
        // infinite loop.
        vector<unsigned int> counts_lo(num_channels, 0);
        // This vector is a rough estimate of what the upper-bound of the
        // pruned_range will be for channel c. This then gets used to call the
        // adjusted_mu() function with reasonable estimates of the numbers of
        // dyes to consider for cross-channel FRET effects. Considerable logic
        // is used for this estimate; the computational burden is low and this
        // improves performance, as the ranges need to be bigger than you would
        // expect.
        vector<unsigned int> counts_hi(num_channels, 0);
        // Repeat until no changes are made to the values. This is set to false
        // on every loop and only flips true when a value is changed.
        bool did_change = true;
        while (did_change) {
            did_change = false;
            for (unsigned int c = 0; c < num_channels; c++) {
                // Convenience variable.
                const ChannelModel& channel_model =
                        *seq_model.channel_models[c];
                // We need to temporarily override the channel c value, so we
                // stash it in another variable then set things how they were.
                // This is how we get the mu value for a count of 1.
                double old_counts_hi_c = counts_hi[c];
                counts_hi[c] = 1;
                double mu1 = channel_model.adjusted_mu(&counts_hi[0]);
                counts_hi[c] = old_counts_hi_c;
                // Need mu value for actual count to estimate sigma (for count).
                double mu = channel_model.adjusted_mu(&counts_hi[0]);
                double s = seq_settings.dist_cutoff * channel_model.sigma(mu);
                // lround rounds to a long (there is no 'round to int' in stl).
                // Also dividing by mu1 is a lazy approximation; due to self-
                // quenching effects, mu grows more slowly as it gets larger.
                // This seems to be a good enough approximation though. Better
                // would be to test different counts until the mu value is as
                // close as possible to the radiometry value.
                unsigned int new_cts = (unsigned int)lround(
                        (radiometry(timestep, c) + s) / mu1);
                // No point in having a value larger than max_num_dyes + 1. We
                // add one because pruned_range max is an exclusive max (<=)
                new_cts = min(max_num_dyes + 1, new_cts);
                // Set and mark changed (or not).
                if (new_cts != counts_hi[c]) {
                    did_change = true;
                    counts_hi[c] = new_cts;
                }
            }
        }
        // Determine pruned_range.
        for (unsigned int c = 0; c < num_channels; c++) {
            // Convenience variable.
            const ChannelModel& channel_model = *seq_model.channel_models[c];

            // Set minimum of pruned range for channel c.
            //
            // min might not get set, if not, the range should be empty, and
            // this is the value we want (pruned_range max is exclusive max).
            unsigned int cmin = max_num_dyes + 1;
            // Stash original value as we will need this to choose cmin for
            // other channels.
            double old_counts_lo_c = counts_lo[c];
            // Increase dye count until the radiometry is in the range (subject
            // to dist_cutoff) of the dye-count. The cut-off is like a z-value
            // in statistics, but is an approximation so that we can more easily
            // account for FRET interactions.
            for (unsigned int d = 0; d < max_num_dyes + 1; d++) {
                counts_lo[c] = d;
                double mu = channel_model.adjusted_mu(&counts_lo[0]);
                double s = seq_settings.dist_cutoff * channel_model.sigma(mu);
                if (mu + s > radiometry(timestep, c)) {
                    cmin = d;
                    break;
                }
            }
            pruned_range.min[1 + c] = cmin;
            // Restore original value of counts_lo[c].
            counts_lo[c] = old_counts_lo_c;

            // Set maximum of pruned range for channel c.
            //
            // max might not get set, if not, this is the value we want
            // (pruned_range max is exclusive max).
            unsigned int cmax = max_num_dyes + 1;
            // Stash original value as we will need this to choose cmax for
            // other channels.
            double old_counts_hi_c = counts_hi[c];
            for (unsigned int d = cmin; d < max_num_dyes + 1; d++) {
                counts_hi[c] = d;
                double mu = channel_model.adjusted_mu(&counts_hi[0]);
                double s = seq_settings.dist_cutoff * channel_model.sigma(mu);
                if (mu - s > radiometry(timestep, c)) {
                    cmax = d;
                    break;
                }
            }
            pruned_range.max[1 + c] = cmax;
            // Restore original value of counts_hi[c].
            counts_hi[c] = old_counts_hi_c;
        }
    }
    // pruned_range saves us the trouble of literally filling in every possible
    // value (with likely performance improvements), but we need to exclude its
    // first element because that refers to the Edman cycle.
    KDRange trange = pruned_range;
    trange.min.erase(trange.min.begin());
    trange.max.erase(trange.max.begin());
    ptsr = new Tensor(trange);
    TensorIterator* it = ptsr->iterator(trange);
    while (!it->done()) {
        *it->get() = 1.0;
        for (unsigned int c = 0; c < num_channels; c++) {
            double observed = radiometry(timestep, c);
            *it->get() *= seq_model.channel_models[c]->pdf(observed, it->loc);
        }
        it->advance();
    }
    delete it;
}

PeptideEmission::PeptideEmission(const PeptideEmission& other)
        : radiometry(other.radiometry),
          timestep(other.timestep),
          pruned_range(other.pruned_range),
          ptsr(other.ptsr),
          i_am_a_copy(true),
          num_channels(other.num_channels),
          max_num_dyes(other.max_num_dyes) {}

PeptideEmission::~PeptideEmission() {
    if (!i_am_a_copy) {
        delete ptsr;
    }
}

void PeptideEmission::prune_forward(KDRange* range, bool* allow_detached) {
    pruned_range = pruned_range.intersect(*range);
    *range = pruned_range;
    this->allow_detached = pruned_range.includes_zero();
    *allow_detached = this->allow_detached;
}

void PeptideEmission::prune_backward(KDRange* range, bool* allow_detached) {
    pruned_range = pruned_range.intersect(*range);
    *range = pruned_range;
    this->allow_detached = pruned_range.includes_zero();
    *allow_detached = this->allow_detached;
}

PeptideStateVector* PeptideEmission::forward_or_backward(
        const PeptideStateVector& input, unsigned int* num_edmans) const {
    vector<unsigned int> zeros(num_channels, 0);
    PeptideStateVector* output = new PeptideStateVector(pruned_range);
    forward_or_backward(input.tensor, &output->tensor);
    forward_or_backward(input.broken_n_tensor, &output->broken_n_tensor);
    if (allow_detached) {
        // This is safe because the allow_detached is only true if pruned_range
        // includes zero.
        double prob = (*ptsr)[&zeros[0]];
        output->p_detached = input.p_detached * prob;
    }
    output->range = pruned_range;
    output->allow_detached = allow_detached;
    return output;
}

void PeptideEmission::forward_or_backward(const Tensor& input,
                                          Tensor* output) const {
    ConstTensorIterator* inputit = input.const_iterator(pruned_range);
    TensorIterator* outputit = output->iterator(pruned_range);
    while (!inputit->done()) {
        // Edman cycle is always the 0th index. We need to snag the rest of the
        // index since the values aren't indexed by Edman cycle.
        double prob = (*ptsr)[&inputit->loc[1]];
        *outputit->get() = *inputit->get() * prob;
        inputit->advance();
        outputit->advance();
    }
    delete inputit;
    delete outputit;
}

PeptideStateVector* PeptideEmission::forward(const PeptideStateVector& input,
                                             unsigned int* num_edmans) const {
    return forward_or_backward(input, num_edmans);
}

PeptideStateVector* PeptideEmission::backward(const PeptideStateVector& input,
                                              unsigned int* num_edmans) const {
    return forward_or_backward(input, num_edmans);
}

void PeptideEmission::improve_fit(const PeptideStateVector& forward_psv,
                                  const PeptideStateVector& backward_psv,
                                  const PeptideStateVector& next_backward_psv,
                                  unsigned int num_edmans,
                                  double probability,
                                  SequencingModelFitter* fitter) const {}

}  // namespace whatprot