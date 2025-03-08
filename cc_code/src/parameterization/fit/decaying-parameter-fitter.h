/******************************************************************************\
* Author: Matthew Beauregard Smith                                             *
* Affiliation: Erisyon, Inc.                                                   *
* Project: Protein Fluorosequencing                                            *
\******************************************************************************/

#ifndef WHATPROT_PARAMETERIZATION_FIT_DECAYING_PARAMETER_FITTER_H
#define WHATPROT_PARAMETERIZATION_FIT_DECAYING_PARAMETER_FITTER_H

// Standard C++ library headers:
#include <vector>

// Local project headers:
#include "parameterization/model/decaying-rate-model.h"
#include "parameterization/settings/fit-settings.h"

namespace whatprot {

class DecayingParameterFitter {
public:
    DecayingParameterFitter();
    DecayingParameterFitter(unsigned int num_timesteps,
                            const DecayingRateModel& prev,
                            const FitSettings& fit_settings);
    void add_timestep(unsigned int t, double x, double n);
    DecayingRateModel get() const;
    DecayingParameterFitter operator+(
            const DecayingParameterFitter& other) const;
    void operator+=(const DecayingParameterFitter& other);
    void operator*=(double weight_adjustment);
    std::vector<double> xvec;
    std::vector<double> nvec;
    DecayingRateModel prev;
    FitSettings fit_settings;
};

}  // namespace whatprot

#endif  // WHATPROT_PARAMETERIZATION_FIT_DECAYING_PARAMETER_FITTER_H
