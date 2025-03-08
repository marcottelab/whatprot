/******************************************************************************\
* Author: Matthew Beauregard Smith                                             *
* Affiliation: The University of Texas at Austin                               *
* Department: Oden Institute and Institute for Cellular and Molecular Biology  *
* PI: Edward Marcotte                                                          *
* Project: Protein Fluorosequencing                                            *
\******************************************************************************/

// Defining symbols from header:
#include "sequencing-model-fitter.h"

// Standard C++ library headers:
#include <utility>

// Local project headers:
#include "parameterization/fit/channel-model-fitter.h"
#include "parameterization/fit/log-normal-distribution-fitter.h"
#include "parameterization/fit/normal-distribution-fitter.h"
#include "parameterization/model/sequencing-model.h"
#include "parameterization/settings/fit-settings.h"

namespace {
using std::move;
}  // namespace

namespace whatprot {

SequencingModelFitter::SequencingModelFitter() {}

SequencingModelFitter::SequencingModelFitter(unsigned int num_timesteps,
                                             unsigned int num_channels,
                                             const SequencingModel& prev,
                                             const FitSettings& fit_settings)
        : p_detach_fit(num_timesteps, prev.p_detach, fit_settings) {
    channel_fits.resize(num_channels);
    for (unsigned int c = 0; c < num_channels; c++) {
        channel_fits[c] = new ChannelModelFitter(*prev.channel_models[c]);
    }
}

SequencingModelFitter::SequencingModelFitter(
        const SequencingModelFitter& other) {
    p_edman_failure_fit = other.p_edman_failure_fit;
    p_detach_fit = other.p_detach_fit;
    p_initial_block_fit = other.p_initial_block_fit;
    p_cyclic_block_fit = other.p_cyclic_block_fit;
    for (unsigned int c = 0; c < other.channel_fits.size(); c++) {
        channel_fits.push_back(new ChannelModelFitter(*other.channel_fits[c]));
    }
}

SequencingModelFitter::SequencingModelFitter(SequencingModelFitter&& other) {
    p_edman_failure_fit = move(other.p_edman_failure_fit);
    p_detach_fit = move(other.p_detach_fit);
    p_initial_block_fit = move(other.p_initial_block_fit);
    p_cyclic_block_fit = move(other.p_cyclic_block_fit);
    channel_fits = move(other.channel_fits);
}

SequencingModelFitter::~SequencingModelFitter() {
    for (ChannelModelFitter* channel_fit : channel_fits) {
        if (channel_fit != NULL) {
            delete channel_fit;
        }
    }
}

SequencingModel SequencingModelFitter::get() const {
    SequencingModel model;
    model.p_edman_failure = p_edman_failure_fit.get();
    model.p_detach = p_detach_fit.get();
    model.p_initial_block = p_initial_block_fit.get();
    model.p_cyclic_block = p_cyclic_block_fit.get();
    for (ChannelModelFitter* channel_fit : channel_fits) {
        model.channel_models.push_back(new ChannelModel(channel_fit->get()));
    }
    return model;
}

SequencingModelFitter SequencingModelFitter::operator+(
        const SequencingModelFitter& other) const {
    SequencingModelFitter result_fitter;
    result_fitter.p_edman_failure_fit =
            p_edman_failure_fit + other.p_edman_failure_fit;
    result_fitter.p_detach_fit = p_detach_fit + other.p_detach_fit;
    result_fitter.p_initial_block_fit =
            p_initial_block_fit + other.p_initial_block_fit;
    result_fitter.p_cyclic_block_fit =
            p_cyclic_block_fit + other.p_cyclic_block_fit;
    for (unsigned int c = 0; c < channel_fits.size(); c++) {
        result_fitter.channel_fits.push_back(new ChannelModelFitter(
                (*channel_fits[c]) + (*other.channel_fits[c])));
    }
    return result_fitter;
}

void SequencingModelFitter::operator+=(const SequencingModelFitter& other) {
    p_edman_failure_fit += other.p_edman_failure_fit;
    p_detach_fit += other.p_detach_fit;
    p_initial_block_fit += other.p_initial_block_fit;
    p_cyclic_block_fit += other.p_cyclic_block_fit;
    for (unsigned int c = 0; c < channel_fits.size(); c++) {
        *channel_fits[c] += *other.channel_fits[c];
    }
}

void SequencingModelFitter::operator*=(double weight_adjustment) {
    p_edman_failure_fit *= weight_adjustment;
    p_detach_fit *= weight_adjustment;
    p_initial_block_fit *= weight_adjustment;
    p_cyclic_block_fit *= weight_adjustment;
    for (unsigned int c = 0; c < channel_fits.size(); c++) {
        *channel_fits[c] *= weight_adjustment;
    }
}

}  // namespace whatprot
