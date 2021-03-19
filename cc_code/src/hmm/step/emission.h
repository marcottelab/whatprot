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
#include "hmm/state-vector/peptide-state-vector.h"
#include "hmm/step/step.h"

namespace whatprot {

class Emission : public Step<PeptideStateVector> {
public:
    Emission(const Radiometry& radiometry,
             int max_num_dyes,
             std::function<double(double, int)> pdf);
    double& prob(int timestep, int channel, int num_dyes);
    double prob(int timestep, int channel, int num_dyes) const;
    virtual void forward(PeptideStateVector* psv) const override;
    virtual void backward(const PeptideStateVector& input,
                          PeptideStateVector* output) const override;
    virtual void improve_fit(const PeptideStateVector& forward_psv,
                             const PeptideStateVector& backward_psv,
                             const PeptideStateVector& next_backward_psv,
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
