/******************************************************************************\
* Author: Matthew Beauregard Smith                                             *
* Affiliation: The University of Texas at Austin                               *
* Department: Oden Institute and Institute for Cellular and Molecular Biology  *
* PI: Edward Marcotte                                                          *
* Project: Protein Fluorosequencing                                            *
\******************************************************************************/

// Defining symbols from header:
#include "channel-model.h"

// Standard C++ library headers:
#include <algorithm>
#include <cmath>
#include <string>

namespace whatprot {

namespace {
using std::abs;
using std::max;
using std::sqrt;
using std::string;
using std::to_string;
}  // namespace

ChannelModel::ChannelModel(unsigned int channel, unsigned int num_channels)
        : channel(channel),
          num_channels(num_channels),
          interactions(num_channels, 1.0) {}

ChannelModel::~ChannelModel() {}

ChannelModel ChannelModel::with_mu_as_one() const {
    ChannelModel x(channel, num_channels);
    x.p_bleach = p_bleach;
    x.p_dud = p_dud;
    x.bg_sig = bg_sig / mu;
    x.mu = 1.0;
    x.sig = sig / mu;
    x.interactions = interactions;
    return x;
}

double ChannelModel::sigma(double amu) const {
    return sqrt(bg_sig * bg_sig + amu * sig * sig);
}

double ChannelModel::distance(const ChannelModel& channel_model) const {
    double dist = 0.0;
    dist = max(dist, abs(p_bleach - channel_model.p_bleach));
    dist = max(dist, abs(p_dud - channel_model.p_dud));
    return dist;
}

double ChannelModel::pdf(double observed, const unsigned int* counts) const {
    return pdf_helper(observed, counts);
}

double ChannelModel::pdf(double observed, const short* counts) const {
    return pdf_helper(observed, counts);
}

string ChannelModel::debug_string() const {
    return "Bleach rate: " + to_string(p_bleach)
           + ", Dud rate: " + to_string(p_dud);
}

}  // namespace whatprot
