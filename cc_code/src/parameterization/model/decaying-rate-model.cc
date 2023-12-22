/******************************************************************************\
* Author: Matthew Beauregard Smith                                             *
* Affiliation: Erisyon, Inc.                                                   *
* Project: Protein Fluorosequencing                                            *
\******************************************************************************/

// Defining symbols from header:
#include "decaying-rate-model.h"

// Standard C++ library headers:
#include <algorithm>
#include <cmath>
#include <utility>

namespace {
using std::abs;
using std::exp;
using std::max;
using std::move;
}  // namespace

namespace whatprot {

DecayingRateModel::DecayingRateModel()
        : base(0), initial(0), initial_decay(0) {}

DecayingRateModel::DecayingRateModel(const DecayingRateModel& other) {
    base = other.base;
    initial = other.initial;
    initial_decay = other.initial_decay;
}

DecayingRateModel& DecayingRateModel::operator=(
        const DecayingRateModel& other) {
    base = other.base;
    initial = other.initial;
    initial_decay = other.initial_decay;
    return *this;
}

DecayingRateModel::DecayingRateModel(DecayingRateModel&& other) {
    base = move(other.base);
    initial = move(other.initial);
    initial_decay = move(other.initial_decay);
}

double DecayingRateModel::operator[](unsigned int i) const {
    // MUST cast to double or garbage output results.
    return base + initial * exp(-(double)i * initial_decay);
}

double DecayingRateModel::distance(const DecayingRateModel& other) const {
    double dist = 0.0;
    dist = max(dist, abs(base - other.base));
    dist = max(dist, abs(initial - other.initial));
    dist = max(dist, abs(initial_decay - other.initial_decay));
    return dist;
}

}  // namespace whatprot
