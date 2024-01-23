/******************************************************************************\
* Author: Matthew Beauregard Smith                                             *
* Affiliation: The University of Texas at Austin                               *
* Department: Oden Institute and Institute for Cellular and Molecular Biology  *
* PI: Edward Marcotte                                                          *
* Project: Protein Fluorosequencing                                            *
\******************************************************************************/

// Defining symbols from header:
#include "channel-model-fitter.h"

// Standard C++ library headers:
#include <utility>

// Local project headers:
#include "parameterization/fit/normal-distribution-fitter.h"
#include "parameterization/model/channel-model.h"

namespace {
using std::move;
}  // namespace

namespace whatprot {

ChannelModelFitter::ChannelModelFitter(const ChannelModel& prev) : prev(prev) {}

ChannelModelFitter::ChannelModelFitter(const ChannelModelFitter& other)
        : prev(other.prev) {
    p_bleach_fit = other.p_bleach_fit;
    p_dud_fit = other.p_dud_fit;
}

ChannelModelFitter::ChannelModelFitter(ChannelModelFitter&& other)
        : prev(other.prev) {
    p_bleach_fit = move(other.p_bleach_fit);
    p_dud_fit = move(other.p_dud_fit);
}

ChannelModel ChannelModelFitter::get() const {
    // Copy variable that for various reasons we do not fit.
    ChannelModel model(prev.channel, prev.num_channels);
    model.interactions = prev.interactions;
    model.bg_sig = prev.bg_sig;
    model.mu = prev.mu;
    model.sig = prev.sig;
    // Then copy fit results.
    model.p_bleach = p_bleach_fit.get();
    model.p_dud = p_dud_fit.get();
    return model;
}

ChannelModelFitter ChannelModelFitter::operator+(
        const ChannelModelFitter& other) const {
    ChannelModelFitter result_fitter(prev);
    result_fitter.p_bleach_fit = p_bleach_fit + other.p_bleach_fit;
    result_fitter.p_dud_fit = p_dud_fit + other.p_dud_fit;
    return result_fitter;
}

void ChannelModelFitter::operator+=(const ChannelModelFitter& other) {
    p_bleach_fit += other.p_bleach_fit;
    p_dud_fit += other.p_dud_fit;
}

void ChannelModelFitter::operator*=(double weight_adjustment) {
    p_bleach_fit *= weight_adjustment;
    p_dud_fit *= weight_adjustment;
}

}  // namespace whatprot
