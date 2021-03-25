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
#include "hmm/fit/error-model-fitter.h"
#include "hmm/state-vector/stuck-dye-state-vector.h"
#include "hmm/step/step.h"

namespace whatprot {

class StuckDyeEmission : public Step<StuckDyeStateVector> {
public:
    StuckDyeEmission(const Radiometry& radiometry,
            int channel,
             std::function<double(double, int)> pdf);
    double& prob(int timestep, int num_dyes);
    double prob(int timestep, int num_dyes) const;
    virtual void forward(StuckDyeStateVector* stsv) const override;
    virtual void backward(const StuckDyeStateVector& input,
                          StuckDyeStateVector* output) const override;
    virtual void improve_fit(const StuckDyeStateVector& forward_stsv,
                             const StuckDyeStateVector& backward_stsv,
                             const StuckDyeStateVector& next_backward_stsv,
                             double probability,
                             ErrorModelFitter* fitter) const override;
    Radiometry radiometry;
    std::vector<double> values;
    int num_timesteps;
    int num_channels;
    int channel;
};

}  // namespace whatprot

#endif  // WHATPROT_HMM_STEP_EMISSION_H
