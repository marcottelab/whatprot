// Author: Matthew Beauregard Smith (UT Austin)
#include "emission.h"

#include <functional>

#include "common/radiometry.h"
#include "tensor/tensor.h"
#include "tensor/tensor_iterator.h"

namespace fluoroseq {

namespace {
using std::function;
}  // namespace

Emission::Emission(const Radiometry& radiometry,
                   int max_num_dyes,
                   function<double (double, int)> pdf)
        : num_timesteps(radiometry.num_timesteps),
          num_channels(radiometry.num_channels),
          max_num_dyes(max_num_dyes) {
    values = new double[num_timesteps * num_channels * max_num_dyes];
    for (int t = 0; t < num_timesteps; t++) {
        for (int c = 0; c < num_channels; c++) {
            for (int d = 0; d < max_num_dyes; d++) {
                prob(t, c, d) = pdf(radiometry(t, c), d);
            }
        }
    }
}

Emission::~Emission() {
    delete[] values;
}

double& Emission::prob(int t, int c, int d) {
    return values[(t * num_channels + c) * max_num_dyes + d];
}

double Emission::prob(int t, int c, int d) const {
    return values[(t * num_channels + c) * max_num_dyes + d];
}

void Emission::operator()(Tensor* tensor, int timestep) const {
    TensorIterator* iterator = tensor->iterator();  // not owned.
    while (!iterator->done()) {
        double product = 1.0;
        for (int c = 0; c < num_channels; c++) {
            product *= prob(timestep, c, iterator->loc[1 + c]);
        }
        *iterator->get() *= product;
        iterator->advance();
    }
}

}  // namespace fluoroseq