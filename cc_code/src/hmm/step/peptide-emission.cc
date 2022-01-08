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
#include <limits>

// Local project headers:
#include "common/radiometry.h"
#include "hmm/state-vector/peptide-state-vector.h"
#include "parameterization/fit/sequencing-model-fitter.h"
#include "parameterization/model/sequencing-model.h"
#include "parameterization/settings/sequencing-settings.h"
#include "tensor/const-tensor-iterator.h"
#include "tensor/tensor-iterator.h"
#include "util/kd-range.h"

namespace whatprot {

namespace {
using std::function;
}  // namespace

PeptideEmission::PeptideEmission(const Radiometry& radiometry,
                                 unsigned int timestep,
                                 int max_num_dyes,
                                 const SequencingModel& seq_model,
                                 const SequencingSettings& seq_settings)
        : radiometry(radiometry),
          timestep(timestep),
          num_channels(radiometry.num_channels),
          max_num_dyes(max_num_dyes) {
    values.resize(num_channels * (max_num_dyes + 1));
    for (unsigned int c = 0; c < num_channels; c++) {
        for (int d = 0; d < (max_num_dyes + 1); d++) {
            prob(c, d) = seq_model.channel_models[c]->pdf(
                    radiometry(timestep, c), d);
        }
    }
    pruned_range.min.resize(1 + num_channels);
    pruned_range.max.resize(1 + num_channels);
    pruned_range.min[0] = 0;
    pruned_range.max[0] = timestep + 1;
    if (seq_settings.dist_cutoff == std::numeric_limits<double>::max()) {
        for (unsigned int c = 0; c < num_channels; c++) {
            pruned_range.min[1 + c] = 0;
            pruned_range.max[1 + c] = std::numeric_limits<unsigned int>::max();
        }
    } else {
        for (unsigned int c = 0; c < num_channels; c++) {
            // min might not get set, if not, the range should be empty, and
            // this is the value we want.
            int cmin = max_num_dyes + 1;
            for (int d = 0; d < max_num_dyes; d++) {
                double s = seq_settings.dist_cutoff
                           * seq_model.channel_models[c]->sigma(d);
                if ((double)d + s > radiometry(timestep, c)) {
                    cmin = d;
                    break;
                }
            }
            pruned_range.min[1 + c] = cmin;
            // max might not get set, if not, this is the value we want.
            int cmax = max_num_dyes + 1;
            for (int d = pruned_range.min[1 + c]; d < max_num_dyes; d++) {
                double s = seq_settings.dist_cutoff
                           * seq_model.channel_models[c]->sigma(d);
                if ((double)d - s > radiometry(timestep, c)) {
                    cmax = d;
                    break;
                }
            }
            pruned_range.max[1 + c] = cmax;
        }
    }
}

double& PeptideEmission::prob(int channel, int num_dyes) {
    return values[channel * (max_num_dyes + 1) + num_dyes];
}

double PeptideEmission::prob(int channel, int num_dyes) const {
    return values[channel * (max_num_dyes + 1) + num_dyes];
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
    PeptideStateVector* output = new PeptideStateVector(pruned_range);
    ConstTensorIterator* inputit = input.tensor.const_iterator(pruned_range);
    TensorIterator* outputit = output->tensor.iterator(pruned_range);
    while (!inputit->done()) {
        double product = 1.0;
        for (unsigned int c = 0; c < num_channels; c++) {
            product *= prob(c, inputit->loc[1 + c]);
        }
        *outputit->get() = *inputit->get() * product;
        inputit->advance();
        outputit->advance();
    }
    delete inputit;
    delete outputit;
    if (allow_detached) {
        double p_detached = input.p_detached;
        for (unsigned int c = 0; c < num_channels; c++) {
            p_detached *= prob(c, 0);
        }
        output->p_detached = p_detached;
    }
    output->range = pruned_range;
    output->allow_detached = allow_detached;
    return output;
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
                                  SequencingModelFitter* fitter) const {
    KDRange range = forward_psv.tensor.range;
    ConstTensorIterator* fit = forward_psv.tensor.const_iterator(range);
    ConstTensorIterator* bit = backward_psv.tensor.const_iterator(range);
    while (fit->index < (num_edmans + 1) * forward_psv.tensor.strides[0]) {
        double p_state = (*fit->get()) * (*bit->get()) / probability;
        for (unsigned int c = 0; c < num_channels; c++) {
            double intensity = radiometry(num_edmans, c);
            int dye_count = fit->loc[1 + c];
            fitter->channel_fits[c]->distribution_fit->add_sample(
                    intensity, dye_count, p_state);
        }
        fit->advance();
        bit->advance();
    }
    delete fit;
    delete bit;
}

}  // namespace whatprot