/******************************************************************************\
* Author: Matthew Beauregard Smith                                             *
* Affiliation: The University of Texas at Austin                               *
* Department: Oden Institute and Institute for Cellular and Molecular Biology  *
* PI: Edward Marcotte                                                          *
* Project: Protein Fluorosequencing                                            *
\******************************************************************************/

#ifndef WHATPROT_HMM_STEP_EMISSION_H
#define WHATPROT_HMM_STEP_EMISSION_H

// Standard C++ library headers:
#include <functional>
#include <vector>

// Local project headers:
#include "common/radiometry.h"
#include "hmm/step/step.h"
#include "tensor/tensor.h"

namespace whatprot {

class Emission : public Step {
public:
    Emission(const Radiometry& radiometry,
             int max_num_dyes,
             std::function<double(double, int)> pdf);
    double& prob(int timestep, int channel, int num_dyes);
    double prob(int timestep, int channel, int num_dyes) const;
    virtual void forward(int* edmans,
                         Tensor* tsr) const override;
    virtual void backward(const Tensor& input,
                          int* edmans,
                          Tensor* output) const override;
    virtual void improve_fit(const Tensor& forward_tensor,
                             const Tensor& backward_tensor,
                             const Tensor& next_backward_tensor,
                             int edmans,
                             double probability,
                             ErrorModelFitter* fitter) const override;
    Radiometry radiometry;
    std::vector<double> values;
    int num_timesteps;
    int num_channels;
    int max_num_dyes;
};

}  // namespace whatprot

#endif  // WHATPROT_HMM_STEP_EMISSION_H
