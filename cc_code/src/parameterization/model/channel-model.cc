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

double ChannelModel::pdf(double observed, int state) const {
    if (state > 0) {
        if (observed == 0.0) {
            return 0.0;
        } else {
            double offset = log(observed) - log((double)state) - mu;
            return (1.0 / (observed * sigma * sqrt(2.0 * PI)))
                   * exp(-(offset * offset) / (2.0 * sigma * sigma));
        }
    } else {
        if (observed == 0.0) {
            return 1.0;
        } else {
            return 0.0;
        }
    }
}

double ChannelModel::relative_distance(
        const ChannelModel& channel_model) const {
    double dist = 0.0;
    dist = max(dist, abs(p_bleach - channel_model.p_bleach) / p_bleach);
    dist = max(dist, abs(p_dud - channel_model.p_dud) / p_dud);
    dist = max(dist, abs(exp(mu) - exp(channel_model.mu)) / exp(mu));
    dist = max(dist, abs(sigma - channel_model.sigma) / sigma);
    dist = max(dist,
               abs(stuck_dye_ratio - channel_model.stuck_dye_ratio)
                       / stuck_dye_ratio);
    dist = max(dist,
               abs(p_stuck_dye_loss - channel_model.p_stuck_dye_loss)
                       / p_stuck_dye_loss);
    return dist;
}

string ChannelModel::debug_string() const {
    return "Bleach rate: " + to_string(p_bleach)
           + ", Dud rate: " + to_string(p_dud)
           + ", exp(mu): " + to_string(exp(mu)) + ", sigma: " + to_string(sigma)
           + ", Stuck dye ratio: " + to_string(stuck_dye_ratio)
           + ", Stuck dye loss rate: " + to_string(p_stuck_dye_loss);
}

}  // namespace whatprot
