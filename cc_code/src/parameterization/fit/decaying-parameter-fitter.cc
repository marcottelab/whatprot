/******************************************************************************\
* Author: Matthew Beauregard Smith                                             *
* Affiliation: Erisyon, Inc.                                                   *
* Project: Protein Fluorosequencing                                            *
\******************************************************************************/

// Defining symbols from header:
#include "decaying-parameter-fitter.h"

// Standard C++ library headers:
#include <vector>

// Local project headers:
#include "parameterization/model/decaying-rate-model.h"
#include "parameterization/settings/fit-settings.h"
#include "util/levmar-helper.h"

namespace {
using std::vector;
}  // namespace

namespace whatprot {

DecayingParameterFitter::DecayingParameterFitter() {}

DecayingParameterFitter::DecayingParameterFitter(
        unsigned int num_timesteps,
        const DecayingRateModel& prev,
        const FitSettings& fit_settings)
        : prev(prev), fit_settings(fit_settings) {
    // One less transition than the number of timesteps - timesteps is the
    // number of emissions.
    xvec.resize(num_timesteps - 1);
    nvec.resize(num_timesteps - 1);
}

void DecayingParameterFitter::add_timestep(unsigned int t, double x, double n) {
    xvec[t] += x;
    nvec[t] += n;
}

DecayingRateModel DecayingParameterFitter::get() const {
    // parameter vector.
    vector<double> params(3);
    params[0] = prev.base;
    params[1] = prev.initial;
    params[2] = prev.initial_decay;
    // y-vector to fit equation to.
    vector<double> yvec(xvec.size());
    for (unsigned int k = 0; k < xvec.size(); k++) {
        yvec[k] = xvec[k] / nvec[k];
    }
    // vector indicating variables to hold.
    vector<bool> hold(3);
    hold[0] = fit_settings.hold_p_detach;
    hold[1] = fit_settings.hold_p_initial_detach;
    hold[2] = fit_settings.hold_p_initial_detach_decay;
    // Run the computation.
    least_squares_fit_of_offset_exponential(yvec, nvec, hold, &params);
    // Extract results and return.
    DecayingRateModel model;
    model.base = params[0];
    model.initial = params[1];
    model.initial_decay = params[2];
    return model;
}

DecayingParameterFitter DecayingParameterFitter::operator+(
        const DecayingParameterFitter& other) const {
    DecayingParameterFitter result_fitter(xvec.size(), prev, fit_settings);
    for (unsigned int i = 0; i < xvec.size(); i++) {
        result_fitter.xvec[i] = xvec[i] + other.xvec[i];
        result_fitter.nvec[i] = nvec[i] + other.nvec[i];
    }
    return result_fitter;
}

void DecayingParameterFitter::operator+=(const DecayingParameterFitter& other) {
    for (unsigned int i = 0; i < xvec.size(); i++) {
        xvec[i] += other.xvec[i];
        nvec[i] += other.nvec[i];
    }
}

void DecayingParameterFitter::operator*=(double weight_adjustment) {
    for (unsigned int i = 0; i < xvec.size(); i++) {
        xvec[i] *= weight_adjustment;
        nvec[i] *= weight_adjustment;
    }
}

}  // namespace whatprot
