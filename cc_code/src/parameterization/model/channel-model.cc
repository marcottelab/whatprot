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

double ChannelModel::pdf(double observed, int state) const {
    double offset = observed - mu * (double)state;
    double s = sigma(state);
    return (1.0 / (s * sqrt(2.0 * PI))) * exp(-offset * offset / (2.0 * s * s));
}

double ChannelModel::sigma(int state) const {
    return sqrt(bg_sig * bg_sig + (double)state * sig * sig);
}

double ChannelModel::relative_distance(
        const ChannelModel& channel_model) const {
    double dist = 0.0;
    dist = max(dist, abs(p_bleach - channel_model.p_bleach) / p_bleach);
    dist = max(dist, abs(p_dud - channel_model.p_dud) / p_dud);
    dist = max(dist, abs(bg_sig - channel_model.bg_sig) / bg_sig);
    dist = max(dist, abs(mu - channel_model.mu) / mu);
    dist = max(dist, abs(sig - channel_model.sig) / sig);
    return dist;
}

string ChannelModel::debug_string() const {
    return "Bleach rate: " + to_string(p_bleach) + ", Dud rate: "
           + to_string(p_dud) + ", mu: " + to_string(mu) + ", sig: "
           + to_string(sig);
}

}  // namespace whatprot
