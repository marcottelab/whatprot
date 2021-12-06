/******************************************************************************\
* Author: Matthew Beauregard Smith                                             *
* Affiliation: The University of Texas at Austin                               *
* Department: Oden Institute and Institute for Cellular and Molecular Biology  *
* PI: Edward Marcotte                                                          *
* Project: Protein Fluorosequencing                                            *
\******************************************************************************/

#ifndef WHATPROT_HMM_STEP_STUCK_DYE_EMISSION_H
#define WHATPROT_HMM_STEP_STUCK_DYE_EMISSION_H

// Standard C++ library headers:
#include <vector>

// Local project headers:
#include "common/radiometry.h"
#include "hmm/state-vector/stuck-dye-state-vector.h"
#include "hmm/step/stuck-dye-step.h"
#include "parameterization/fit/sequencing-model-fitter.h"
#include "parameterization/model/sequencing-model.h"

namespace whatprot {

class StuckDyeEmission : public StuckDyeStep {
public:
    StuckDyeEmission(const Radiometry& radiometry,
                     int channel,
                     const SequencingModel& seq_model);
    double& prob(int timestep, int channel, int num_dyes);
    double prob(int timestep, int channel, int num_dyes) const;
    virtual void forward(unsigned int* num_edmans,
                         StuckDyeStateVector* sdsv) const override;
    virtual void backward(const StuckDyeStateVector& input,
                          unsigned int* num_edmans,
                          StuckDyeStateVector* output) const override;
    virtual void improve_fit(const StuckDyeStateVector& forward_sdsv,
                             const StuckDyeStateVector& backward_sdsv,
                             const StuckDyeStateVector& next_backward_sdsv,
                             unsigned int num_edmans,
                             double probability,
                             SequencingModelFitter* fitter) const override;
    Radiometry radiometry;
    std::vector<double> values;
    unsigned int num_timesteps;
    unsigned int num_channels;
    unsigned int channel;
};

}  // namespace whatprot

#endif  // WHATPROT_HMM_STEP_STUCK_DYE_EMISSION_H
