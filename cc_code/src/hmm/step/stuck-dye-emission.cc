/******************************************************************************\
* Author: Matthew Beauregard Smith                                             *
* Affiliation: The University of Texas at Austin                               *
* Department: Oden Institute and Institute for Cellular and Molecular Biology  *
* PI: Edward Marcotte                                                          *
* Project: Protein Fluorosequencing                                            *
\******************************************************************************/

// Defining symbols from header:
#include "stuck-dye-emission.h"

// Standard C++ library headers:
#include <functional>

// Local project headers:
#include "common/radiometry.h"
#include "hmm/fit/error-model-fitter.h"
#include "hmm/state-vector/peptide-state-vector.h"

namespace whatprot {

namespace {
using std::function;
}  // namespace

StuckDyeEmission::StuckDyeEmission(const Radiometry& radiometry,
                                   int channel,
                                   function<double(double, int)> pdf)
        : radiometry(radiometry),
          num_timesteps(radiometry.num_timesteps),
          num_channels(radiometry.num_channels),
          channel(channel) {
    values.resize(num_timesteps * num_channels * 2);
    for (int t = 0; t < num_timesteps; t++) {
        for (int c = 0; c < num_channels; c++) {
            for (int d = 0; d < 2; d++) {
                prob(t, c, d) = pdf(radiometry(t, c), d);
            }
        }
    }
}

double& StuckDyeEmission::prob(int t, int c, int d) {
    return values[(t * num_channels + c) * 2 + d];
}

double StuckDyeEmission::prob(int t, int c, int d) const {
    return values[(t * num_channels + c) * 2 + d];
}

void StuckDyeEmission::forward(int* num_edmans,
                               StuckDyeStateVector* sdsv) const {
    double product = 1.0;
    for (int c = 0; c < num_channels; c++) {
        if (c == channel) {
            continue;
        }
        product *= prob(*num_edmans, c, 0);
    }
    sdsv->dye *= product * prob(*num_edmans, channel, 1);
    sdsv->no_dye *= product * prob(*num_edmans, channel, 0);
}

void StuckDyeEmission::backward(const StuckDyeStateVector& input,
                                int* num_edmans,
                                StuckDyeStateVector* output) const {
    double product = 1.0;
    for (int c = 0; c < num_channels; c++) {
        if (c == channel) {
            continue;
        }
        product *= prob((*num_edmans), c, 0);
    }
    output->dye *= input.dye * product * prob((*num_edmans), channel, 1);
    output->no_dye *= input.no_dye * product * prob((*num_edmans), channel, 0);
}

void StuckDyeEmission::improve_fit(
        const StuckDyeStateVector& forward_sdsv,
        const StuckDyeStateVector& backward_sdsv,
        const StuckDyeStateVector& next_backward_sdsv,
        int num_edmans,
        double probability,
        ErrorModelFitter* fitter) const {
    double intensity = radiometry(num_edmans, channel);
    double p_no_dye = forward_sdsv.no_dye * backward_sdsv.no_dye;
    fitter->distribution_fit->add_sample(intensity, 0, p_no_dye);
    double p_dye = forward_sdsv.dye * backward_sdsv.dye;
    fitter->distribution_fit->add_sample(intensity, 1, p_dye);
    double p_total = p_no_dye + p_dye;
    for (int c = 0; c < num_channels; c++) {
        if (c == channel) {
            continue;
        }
        intensity = radiometry(num_edmans, c);
        fitter->distribution_fit->add_sample(intensity, 0, p_total);
    }
}

}  // namespace whatprot