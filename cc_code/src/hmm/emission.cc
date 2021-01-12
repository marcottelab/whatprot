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

namespace fluoroseq {

namespace {
using std::function;
}  // namespace

Emission::Emission(const Radiometry& radiometry,
                   int max_num_dyes,
                   function<double(double, int)> pdf)
        : num_timesteps(radiometry.num_timesteps),
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

void Emission::forward(const Tensor& input, int timestep, Tensor* output) const {
    ConstTensorIterator* inputit = input.const_iterator();
    TensorIterator* outputit = output->iterator();
    while (inputit->index < (timestep + 1) * input.strides[0]) {
        double product = 1.0;
        for (int c = 0; c < num_channels; c++) {
            product *= prob(timestep, c, inputit->loc[1 + c]);
        }
        *outputit->get() = inputit->get() * product;
        inputit->advance();
        outputit->advance();
    }
    delete inputit;
    delete outputit;
}

}  // namespace fluoroseq