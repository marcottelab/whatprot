/******************************************************************************\
* Author: Matthew Beauregard Smith                                             *
* Affiliation: The University of Texas at Austin                               *
* Department: Oden Institute and Institute for Cellular and Molecular Biology  *
* PI: Edward Marcotte                                                          *
* Project: Protein Fluorosequencing                                            *
\******************************************************************************/

// Defining symbols from header:
#include "peptide-emission.h"

// Local project headers:
#include "common/radiometry.h"
#include "hmm/state-vector/peptide-state-vector.h"
#include "parameterization/fit/sequencing-model-fitter.h"
#include "parameterization/model/sequencing-model.h"
#include "tensor/const-tensor-iterator.h"
#include "tensor/tensor-iterator.h"
#include "util/kd-box-range.h"

namespace whatprot {

namespace {
using std::function;
}  // namespace

PeptideEmission::PeptideEmission(const Radiometry& radiometry,
                                 unsigned int timestep,
                                 int max_num_dyes,
                                 const SequencingModel& seq_model)
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
}

double& PeptideEmission::prob(int channel, int num_dyes) {
    return values[channel * (max_num_dyes + 1) + num_dyes];
}

double PeptideEmission::prob(int channel, int num_dyes) const {
    return values[channel * (max_num_dyes + 1) + num_dyes];
}

void PeptideEmission::forward(unsigned int* num_edmans,
                              PeptideStateVector* psv) const {
    KDBoxRange range;
    range.min.resize(1 + num_channels);
    range.max.resize(1 + num_channels);
    for (unsigned int o = 0; o < 1 + num_channels; o++) {
        range.min[o] = 0;
        range.max[o] = psv->tensor.shape[o];
    }
    TensorIterator* it = psv->tensor.iterator(range);
    while (it->index < (*num_edmans + 1) * psv->tensor.strides[0]) {
        double product = 1.0;
        for (unsigned int c = 0; c < num_channels; c++) {
            product *= prob(c, it->loc[1 + c]);
        }
        *it->get() = *it->get() * product;
        it->advance();
    }
    delete it;
}

void PeptideEmission::backward(const PeptideStateVector& input,
                               unsigned int* num_edmans,
                               PeptideStateVector* output) const {
    KDBoxRange range;
    range.min.resize(1 + num_channels);
    range.max.resize(1 + num_channels);
    for (unsigned int o = 0; o < 1 + num_channels; o++) {
        range.min[o] = 0;
        range.max[o] = input.tensor.shape[o];
    }
    ConstTensorIterator* inputit = input.tensor.const_iterator(range);
    TensorIterator* outputit = output->tensor.iterator(range);
    while (inputit->index < (*num_edmans + 1) * input.tensor.strides[0]) {
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
}

void PeptideEmission::improve_fit(const PeptideStateVector& forward_psv,
                                  const PeptideStateVector& backward_psv,
                                  const PeptideStateVector& next_backward_psv,
                                  unsigned int num_edmans,
                                  double probability,
                                  SequencingModelFitter* fitter) const {
    KDBoxRange range;
    range.min.resize(1 + num_channels);
    range.max.resize(1 + num_channels);
    for (unsigned int o = 0; o < 1 + num_channels; o++) {
        range.min[o] = 0;
        range.max[o] = forward_psv.tensor.shape[o];
    }
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