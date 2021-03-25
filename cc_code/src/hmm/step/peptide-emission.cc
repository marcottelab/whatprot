/******************************************************************************\
* Author: Matthew Beauregard Smith                                             *
* Affiliation: The University of Texas at Austin                               *
* Department: Oden Institute and Institute for Cellular and Molecular Biology  *
* PI: Edward Marcotte                                                          *
* Project: Protein Fluorosequencing                                            *
\******************************************************************************/

// Defining symbols from header:
#include "emission.h"

// Standard C++ library headers:
#include <functional>

// Local project headers:
#include "common/radiometry.h"
#include "hmm/fit/error-model-fitter.h"
#include "hmm/state-vector/peptide-state-vector.h"
#include "tensor/tensor-iterator.h"

namespace whatprot {

namespace {
using std::function;
}  // namespace

PeptideEmission::PeptideEmission(const Radiometry& radiometry,
                   int max_num_dyes,
                   function<double(double, int)> pdf)
        : radiometry(radiometry),
          num_timesteps(radiometry.num_timesteps),
          num_channels(radiometry.num_channels),
          max_num_dyes(max_num_dyes) {
    values.resize(num_timesteps * num_channels * (max_num_dyes + 1));
    for (int t = 0; t < num_timesteps; t++) {
        for (int c = 0; c < num_channels; c++) {
            for (int d = 0; d < (max_num_dyes + 1); d++) {
                prob(t, c, d) = pdf(radiometry(t, c), d);
            }
        }
    }
}

double& PeptideEmission::prob(int t, int c, int d) {
    return values[(t * num_channels + c) * (max_num_dyes + 1) + d];
}

double PeptideEmission::prob(int t, int c, int d) const {
    return values[(t * num_channels + c) * (max_num_dyes + 1) + d];
}

void PeptideEmission::forward(PeptideStateVector* psv) const {
    TensorIterator* it = psv->tensor.iterator();
    while (it->index < (psv->num_edmans + 1) * psv->tensor.strides[0]) {
        double product = 1.0;
        for (int c = 0; c < num_channels; c++) {
            product *= prob(psv->num_edmans, c, it->loc[1 + c]);
        }
        *it->get() = *it->get() * product;
        it->advance();
    }
    delete it;
}

void PeptideEmission::backward(const PeptideStateVector& input,
                        PeptideStateVector* output) const {
    ConstTensorIterator* inputit = input.tensor.const_iterator();
    TensorIterator* outputit = output->tensor.iterator();
    while (inputit->index < (input.num_edmans + 1) * input.tensor.strides[0]) {
        double product = 1.0;
        for (int c = 0; c < num_channels; c++) {
            product *= prob(input.num_edmans, c, inputit->loc[1 + c]);
        }
        *outputit->get() = inputit->get() * product;
        inputit->advance();
        outputit->advance();
    }
    delete inputit;
    delete outputit;
    output->num_edmans = input.num_edmans;
}

void PeptideEmission::improve_fit(const PeptideStateVector& forward_psv,
                           const PeptideStateVector& backward_psv,
                           const PeptideStateVector& next_backward_psv,
                           double probability,
                           ErrorModelFitter* fitter) const {
    ConstTensorIterator* fit = forward_psv.tensor.const_iterator();
    ConstTensorIterator* bit = backward_psv.tensor.const_iterator();
    while (fit->index
           < (forward_psv.num_edmans + 1) * forward_psv.tensor.strides[0]) {
        double p_state = fit->get() * bit->get() / probability;
        for (int c = 0; c < num_channels; c++) {
            int t = fit->loc[0];
            double intensity = radiometry(forward_psv.num_edmans, c);
            int dye_count = fit->loc[1 + c];
            fitter->distribution_fit->add_sample(intensity, dye_count, p_state);
        }
        fit->advance();
        bit->advance();
    }
    delete fit;
    delete bit;
}

}  // namespace whatprot