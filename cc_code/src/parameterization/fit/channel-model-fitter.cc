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

ChannelModelFitter::ChannelModelFitter() {
    distribution_fit = new NormalDistributionFitter();
}

ChannelModelFitter::ChannelModelFitter(const ChannelModelFitter& other) {
    p_initial_bleach_fit = other.p_initial_bleach_fit;
    p_cyclic_bleach_fit = other.p_cyclic_bleach_fit;
    p_dud_fit = other.p_dud_fit;
    distribution_fit = new NormalDistributionFitter(*other.distribution_fit);
}

ChannelModelFitter::ChannelModelFitter(ChannelModelFitter&& other) {
    p_initial_bleach_fit = move(other.p_initial_bleach_fit);
    p_cyclic_bleach_fit = move(other.p_cyclic_bleach_fit);
    p_dud_fit = move(other.p_dud_fit);
    distribution_fit = other.distribution_fit;
    other.distribution_fit = NULL;
}

ChannelModelFitter::~ChannelModelFitter() {
    if (distribution_fit != NULL) {
        delete distribution_fit;
    }
}

ChannelModel ChannelModelFitter::get() const {
    ChannelModel model;
    model.p_initial_bleach = p_initial_bleach_fit.get();
    model.p_cyclic_bleach = p_cyclic_bleach_fit.get();
    model.p_dud = p_dud_fit.get();
    model.bg_sig = distribution_fit->get_bg_sig();
    model.mu = distribution_fit->get_mu();
    model.sig = distribution_fit->get_sig();
    return model;
}

ChannelModelFitter ChannelModelFitter::operator+(
        const ChannelModelFitter& other) const {
    ChannelModelFitter result_fitter;
    result_fitter.p_initial_bleach_fit =
            p_initial_bleach_fit + other.p_initial_bleach_fit;
    result_fitter.p_cyclic_bleach_fit =
            p_cyclic_bleach_fit + other.p_cyclic_bleach_fit;
    result_fitter.p_dud_fit = p_dud_fit + other.p_dud_fit;
    *result_fitter.distribution_fit =
            *distribution_fit + *other.distribution_fit;
    return result_fitter;
}

void ChannelModelFitter::operator+=(const ChannelModelFitter& other) {
    p_initial_bleach_fit += other.p_initial_bleach_fit;
    p_cyclic_bleach_fit += other.p_cyclic_bleach_fit;
    p_dud_fit += other.p_dud_fit;
    *distribution_fit += *other.distribution_fit;
}

void ChannelModelFitter::operator*=(double weight_adjustment) {
    p_initial_bleach_fit *= weight_adjustment;
    p_cyclic_bleach_fit *= weight_adjustment;
    p_dud_fit *= weight_adjustment;
    *distribution_fit *= weight_adjustment;
}

}  // namespace whatprot
