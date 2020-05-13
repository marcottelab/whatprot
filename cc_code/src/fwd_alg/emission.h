// Author: Matthew Beauregard Smith (UT Austin)
#ifndef FLUOROSEQ_FWD_ALG_EMISSION
#define FLUOROSEQ_FWD_ALG_EMISSION

#include <functional>

#include "common/radiometry.h"
#include "tensor/tensor.h"

namespace fluoroseq {

namespace {
using std::function;
}  // namespace

class Emission {
public:
    Emission(const Radiometry& radiometry,
             int max_num_dyes,
             function<double (double, int)> pdf,
             int max_edman_failures);
    ~Emission();
    double& prob(int timestep, int channel, int num_dyes);
    double prob(int timestep, int channel, int num_dyes) const;
    void operator()(Tensor* tensor, int timestep) const;

    double* values;
    int num_timesteps;
    int num_channels;
    int max_num_dyes;
    int max_edman_failures;
};

}  // namespace fluoroseq

#endif  // FLUOROSEQ_FWD_ALG_EMISSION
