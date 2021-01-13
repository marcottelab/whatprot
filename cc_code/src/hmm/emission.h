/******************************************************************************\
* Author: Matthew Beauregard Smith                                             *
* Affiliation: The University of Texas at Austin                               *
* Department: Oden Institute and Institute for Cellular and Molecular Biology  *
* PI: Edward Marcotte                                                          *
* Project: Protein Fluorosequencing                                            *
\******************************************************************************/

#ifndef FLUOROSEQ_HMM_EMISSION
#define FLUOROSEQ_HMM_EMISSION

// Standard C++ library headers:
#include <functional>
#include <vector>

// Local project headers:
#include "common/radiometry.h"
#include "tensor/tensor.h"

namespace fluoroseq {

class Emission {
public:
    Emission(const Radiometry& radiometry,
             int max_num_dyes,
             std::function<double(double, int)> pdf);
    double& prob(int timestep, int channel, int num_dyes);
    double prob(int timestep, int channel, int num_dyes) const;
    void forward(const Tensor& input, int timestep, Tensor* output) const;
    void backward(const Tensor& input, int timestep, Tensor* output) const;

    std::vector<double> values;
    int num_timesteps;
    int num_channels;
    int max_num_dyes;
};

}  // namespace fluoroseq

#endif  // FLUOROSEQ_HMM_EMISSION
