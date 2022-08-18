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
#include "hmm/state-vector/peptide-state-vector.h"
#include "parameterization/fit/sequencing-model-fitter.h"
#include "parameterization/model/sequencing-model.h"

namespace whatprot {

namespace {
using std::function;
}  // namespace

StuckDyeEmission::StuckDyeEmission(const Radiometry& radiometry,
                                   int channel,
                                   const SequencingModel& seq_model,
                    const SequencingSettings& seq_settings)
        : radiometry(radiometry),
          num_timesteps(radiometry.num_timesteps),
          num_channels(radiometry.num_channels),
          channel(channel) {
    double s0 = seq_settings.dist_cutoff
                           * seq_model.channel_models[channel]->sigma(0);
    double s1 = seq_settings.dist_cutoff
                           * seq_model.channel_models[channel]->sigma(1);
    double m1 = seq_model.channel_models[channel]->mu;
    values.resize(num_timesteps * num_channels * 2);
    for (unsigned int t = 0; t < num_timesteps; t++) {
        for (unsigned int c = 0; c < num_channels; c++) {
            double x = radiometry(t, c);
            if (-s0 < x && x < s0) {
            prob(t, c, 0) = seq_model.channel_models[c]->pdf(x, 0);
            } else {
                prob(t, c, 0) = 0.0;
            }
            if (m1 - s1 < x && x < m1 + s1) {
                prob(t, c, 1) = seq_model.channel_models[c]->pdf(x, 0);
            }else {
                prob(t, c, 1) = 0.0;
            }
            // for (int d = 0; d < 2; d++) {
            //     prob(t, c, d) =
            //             seq_model.channel_models[c]->pdf(radiometry(t, c), d);
            // }
        }
    }
}

double& StuckDyeEmission::prob(int t, int c, int d) {
    return values[(t * num_channels + c) * 2 + d];
}

double StuckDyeEmission::prob(int t, int c, int d) const {
    return values[(t * num_channels + c) * 2 + d];
}

StuckDyeStateVector* StuckDyeEmission::forward(const StuckDyeStateVector& input,
                                               unsigned int* num_edmans) const {
    double product = 1.0;
    for (unsigned int c = 0; c < num_channels; c++) {
        if (c == channel) {
            continue;
        }
        product *= prob(*num_edmans, c, 0);
    }
    StuckDyeStateVector* output = new StuckDyeStateVector();
    output->dye = input.dye * product * prob(*num_edmans, channel, 1);
    output->no_dye = input.no_dye * product * prob(*num_edmans, channel, 0);
    return output;
}

StuckDyeStateVector* StuckDyeEmission::backward(
        const StuckDyeStateVector& input, unsigned int* num_edmans) const {
    return forward(input, num_edmans);
}

void StuckDyeEmission::improve_fit(
        const StuckDyeStateVector& forward_sdsv,
        const StuckDyeStateVector& backward_sdsv,
        const StuckDyeStateVector& next_backward_sdsv,
        unsigned int num_edmans,
        double probability,
        SequencingModelFitter* fitter) const {
    double intensity = radiometry(num_edmans, channel);
    double p_no_dye = forward_sdsv.no_dye * backward_sdsv.no_dye / probability;
    fitter->channel_fits[channel]->distribution_fit->add_sample(
            intensity, 0, p_no_dye);
    double p_dye = forward_sdsv.dye * backward_sdsv.dye / probability;
    fitter->channel_fits[channel]->distribution_fit->add_sample(
            intensity, 1, p_dye);
    double p_total = p_no_dye + p_dye;
    for (unsigned int c = 0; c < num_channels; c++) {
        if (c == channel) {
            continue;
        }
        intensity = radiometry(num_edmans, c);
        fitter->channel_fits[c]->distribution_fit->add_sample(
                intensity, 0, p_total);
    }
}

}  // namespace whatprot