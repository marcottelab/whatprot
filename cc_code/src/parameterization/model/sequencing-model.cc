/******************************************************************************\
* Author: Matthew Beauregard Smith                                             *
* Affiliation: The University of Texas at Austin                               *
* Department: Oden Institute and Institute for Cellular and Molecular Biology  *
* PI: Edward Marcotte                                                          *
* Project: Protein Fluorosequencing                                            *
\******************************************************************************/

// Defining symbols from header:
#include "sequencing-model.h"

// Standard C++ library headers:
#include <algorithm>
#include <cmath>
#include <string>

// Local project headers:
#include "parameterization/model/channel-model.h"

namespace {
using std::abs;
using std::exp;
using std::max;
using std::move;
using std::string;
using std::to_string;
}  // namespace

namespace whatprot {

SequencingModel::SequencingModel() {}

SequencingModel::SequencingModel(const SequencingModel& other) {
    p_edman_failure = other.p_edman_failure;
    p_detach = other.p_detach;
    for (unsigned int c = 0; c < other.channel_models.size(); c++) {
        channel_models.push_back(new ChannelModel(*other.channel_models[c]));
    }
}

SequencingModel& SequencingModel::operator=(const SequencingModel& other) {
    p_edman_failure = other.p_edman_failure;
    p_detach = other.p_detach;
    for (unsigned int c = 0; c < other.channel_models.size(); c++) {
        channel_models.push_back(new ChannelModel(*other.channel_models[c]));
    }
    return *this;
}

SequencingModel::SequencingModel(SequencingModel&& other) {
    p_edman_failure = move(other.p_edman_failure);
    p_detach = move(other.p_detach);
    channel_models = move(other.channel_models);
}

SequencingModel::~SequencingModel() {
    for (ChannelModel* channel_model : channel_models) {
        if (channel_model != NULL) {
            delete channel_model;
        }
    }
}

double SequencingModel::relative_distance(
        const SequencingModel& sequencing_model) const {
    double dist = 0.0;
    dist = max(dist,
               abs(p_edman_failure - sequencing_model.p_edman_failure)
                       / p_edman_failure);
    dist = max(dist, abs(p_detach - sequencing_model.p_detach) / p_detach);
    for (unsigned int i = 0; i < channel_models.size(); i++) {
        dist = max(dist,
                   channel_models[i]->relative_distance(
                           *sequencing_model.channel_models[i]));
    }
    return dist;
}

string SequencingModel::debug_string() const {
    string s = "";
    s += "Edman failure rate: " + to_string(p_edman_failure) + ", ";
    s += "Detach rate: " + to_string(p_detach);
    for (unsigned int i = 0; i < channel_models.size(); i++) {
        s += ", ";
        s += "Channel " + to_string(i) + " info: (";
        s += channel_models[i]->debug_string() + ")";
    }
    return s;
}

}  // namespace whatprot
