// Author: Matthew Beauregard Smith (UT Austin)
#ifndef FLUOROSEQ_FWD_ALG_EMISSION
#define FLUOROSEQ_FWD_ALG_EMISSION

#include <functional>
#include <vector>

#include "common/radiometry.h"
#include "tensor/tensor.h"

namespace fluoroseq {

class Emission {
public:
    Emission(const Radiometry& radiometry,
             int max_num_dyes,
             std::function<double (double, int)> pdf);
    double& prob(int timestep, int channel, int num_dyes);
    double prob(int timestep, int channel, int num_dyes) const;
    void operator()(Tensor* tensor, int timestep) const;

    std::vector<double> values;
    int num_timesteps;
    int num_channels;
    int max_num_dyes;
};

}  // namespace fluoroseq

#endif  // FLUOROSEQ_FWD_ALG_EMISSION
