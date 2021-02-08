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
#include "tensor/tensor-iterator.h"
#include "tensor/tensor.h"

namespace whatprot {

namespace {
using std::function;
}  // namespace

Emission::Emission(const Radiometry& radiometry,
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

double& Emission::prob(int t, int c, int d) {
    return values[(t * num_channels + c) * (max_num_dyes + 1) + d];
}

double Emission::prob(int t, int c, int d) const {
    return values[(t * num_channels + c) * (max_num_dyes + 1) + d];
}

void Emission::forward(const Tensor& input, int* edmans, Tensor* output) const {
    ConstTensorIterator* inputit = input.const_iterator();
    TensorIterator* outputit = output->iterator();
    while (inputit->index < (*edmans + 1) * input.strides[0]) {
        double product = 1.0;
        for (int c = 0; c < num_channels; c++) {
            product *= prob(*edmans, c, inputit->loc[1 + c]);
        }
        *outputit->get() = inputit->get() * product;
        inputit->advance();
        outputit->advance();
    }
    delete inputit;
    delete outputit;
}

void Emission::backward(const Tensor& input,
                        int* edmans,
                        Tensor* output) const {
    forward(input, edmans, output);
}

void Emission::improve_fit(const Tensor& forward_tensor,
                           const Tensor& backward_tensor,
                           const Tensor& next_backward_tensor,
                           int edmans,
                           double probability,
                           ErrorModelFitter* fitter) const {
    ConstTensorIterator* fit = forward_tensor.const_iterator();
    ConstTensorIterator* bit = backward_tensor.const_iterator();
    while (fit->index < (edmans + 1) * forward_tensor.strides[0]) {
        double p_state = fit->get() * bit->get() / probability;
        for (int c = 0; c < num_channels; c++) {
            int t = fit->loc[0];
            double intensity = radiometry(t, c);
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