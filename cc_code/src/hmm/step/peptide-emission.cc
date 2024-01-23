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
        for (unsigned int c = 0; c < num_channels; c++) {
            // We use this a lot in this scope. Nice as a convenience variable.
            const ChannelModel& channel_model = *seq_model.channel_models[c];
            // Need this to collect adjusted_mu (which is in turn to get sigma).
            vector<unsigned int> counts(num_channels, 0);
            // min might not get set, if not, the range should be empty, and
            // this is the value we want.
            unsigned int cmin = max_num_dyes + 1;
            for (unsigned int d = 0; d < max_num_dyes; d++) {
                counts[c] = d;
                double mu = channel_model.adjusted_mu(&counts[0]);
                double s = seq_settings.dist_cutoff * channel_model.sigma(mu);
                if ((double)d * mu + s > radiometry(timestep, c)) {
                    cmin = d;
                    break;
                }
            }
            pruned_range.min[1 + c] = cmin;
            // Need this at non-zero value in order to test effect on mu of
            // cross-channel interactions.
            counts[c] = 1;
            // To set maximum end of the pruned ranges we need the lower bound
            // of distributions. This is variable due to cross-dye interactions.
            // Extreme scenarios exist in some cases with max_num_dyes larger
            // than 100. To avoid overly impacting classification against other
            // peptides, we use cross_dye_maximum so that we can avoid
            // overpruning.
            for (unsigned int i = 0; i < seq_settings.cross_dye_maximum; i++) {
                unsigned int cross_channel = 0;
                bool has_interaction = false;
                double worst_mu = channel_model.adjusted_mu(&counts[0]);
                for (unsigned int cc = 0; cc < num_channels; cc++) {
                    // Skip self. Self quenching effect is handled gracefully by
                    // the adjusted mu calculation without creating the same
                    // kind of problems too many fluorophores does.
                    if (cc == c) {
                        continue;
                    }
                    // Never put more dyes on a channel than the max_num_dyes.
                    if (counts[cc] == max_num_dyes) {
                        continue;
                    }
                    // Easier to use existing variable than create new one, so
                    // we increase channel cc, see what value it gives, then
                    // return it to its original value.
                    counts[cc] += 1;
                    double mu = channel_model.adjusted_mu(&counts[0]);
                    counts[cc] -= 1;
                    if (mu < worst_mu) {
                        cross_channel = cc;
                        worst_mu = mu;
                        has_interaction = true;
                    }
                }
                // If we didn't find an interaction each loop will be the same
                // and we never will.
                if (!has_interaction) {
                    break;
                }
                // Finally we can update the value.
                counts[cross_channel] += 1;
            }
            // max might not get set, if not, this is the value we want.
            unsigned int cmax = max_num_dyes + 1;
            for (unsigned int d = cmin; d < max_num_dyes; d++) {
                counts[c] = d;
                double mu = channel_model.adjusted_mu(&counts[0]);
                double s = seq_settings.dist_cutoff * channel_model.sigma(mu);
                if ((double)d - s > radiometry(timestep, c)) {
                    cmax = d;
                    break;
                }
            }
            pruned_range.max[1 + c] = cmax;
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
