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
#include <functional>
#include <string>

namespace whatprot {

namespace {
using std::abs;
using std::exp;
using std::function;
using std::log;
using std::max;
using std::sqrt;
using std::string;
using std::to_string;
double PI = 3.141592653589793238;
}  // namespace

ChannelModel::~ChannelModel() {}

ChannelModel ChannelModel::with_mu_as_one() const {
    ChannelModel x;
    x.p_initial_bleach = p_initial_bleach;
    x.p_cyclic_bleach = p_cyclic_bleach;
    x.p_dud = p_dud;
    x.bg_sig = bg_sig / mu;
    x.mu = 1.0;
    x.sig = sig / mu;
    return x;
}

double ChannelModel::pdf(double observed, int state) const {
    double offset = observed - mu * (double)state;
    double s = sigma(state);
    return (1.0 / (s * sqrt(2.0 * PI))) * exp(-offset * offset / (2.0 * s * s));
}

double ChannelModel::sigma(int state) const {
    return sqrt(bg_sig * bg_sig + (double)state * sig * sig);
}

double ChannelModel::distance(const ChannelModel& channel_model) const {
    double dist = 0.0;
    dist = max(dist, abs(p_initial_bleach - channel_model.p_initial_bleach));
    dist = max(dist, abs(p_cyclic_bleach - channel_model.p_cyclic_bleach));
    dist = max(dist, abs(p_dud - channel_model.p_dud));
    dist = max(dist, abs(bg_sig - channel_model.bg_sig));
    dist = max(dist, abs(mu - channel_model.mu));
    dist = max(dist, abs(sig - channel_model.sig));
    return dist;
}

string ChannelModel::debug_string() const {
    return "Initial bleach rate: " + to_string(p_initial_bleach)
           + ", Cyclic bleach rate: " + to_string(p_cyclic_bleach)
           + ", Dud rate: " + to_string(p_dud);
}

}  // namespace whatprot
