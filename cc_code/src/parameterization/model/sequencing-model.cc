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
#include <fstream>
#include <string>

// External headers:
#include "json.hpp"

// Local project headers:
#include "parameterization/model/channel-model.h"

namespace {
using json = nlohmann::json;
using std::abs;
using std::exp;
using std::ifstream;
using std::max;
using std::move;
using std::string;
using std::to_string;
}  // namespace

namespace whatprot {

SequencingModel::SequencingModel() {}

SequencingModel::SequencingModel(const string& seq_model_filename) {
    ifstream f(seq_model_filename);
    json data = json::parse(f);
    p_edman_failure = data["p_edman_failure"].get<double>();
    p_detach = data["p_detach"].get<double>();
    p_initial_break_n = data["p_initial_break_n"].get<double>();
    p_cyclic_break_n = data["p_cyclic_break_n"].get<double>();
    for (auto& channel_data : data["channel_models"]) {
        channel_models.push_back(new ChannelModel());
        channel_models.back()->p_bleach =
                channel_data["p_bleach"].get<double>();
        channel_models.back()->p_dud = channel_data["p_dud"].get<double>();
        channel_models.back()->bg_sig = channel_data["bg_sig"].get<double>();
        channel_models.back()->mu = channel_data["mu"].get<double>();
        channel_models.back()->sig = channel_data["sig"].get<double>();
    }
}

SequencingModel::SequencingModel(const SequencingModel& other) {
    p_edman_failure = other.p_edman_failure;
    p_detach = other.p_detach;
    p_initial_break_n = other.p_initial_break_n;
    p_cyclic_break_n = other.p_cyclic_break_n;
    for (unsigned int c = 0; c < other.channel_models.size(); c++) {
        channel_models.push_back(new ChannelModel(*other.channel_models[c]));
    }
}

SequencingModel& SequencingModel::operator=(const SequencingModel& other) {
    p_edman_failure = other.p_edman_failure;
    p_detach = other.p_detach;
    p_initial_break_n = other.p_initial_break_n;
    p_cyclic_break_n = other.p_cyclic_break_n;
    // This function is not necessarily used as a constructor. It is very
    // important to clear contents of channel_models before filling it.
    for (ChannelModel* channel_model : channel_models) {
        if (channel_model != NULL) {
            delete channel_model;
        }
    }
    channel_models.resize(0);
    // Now we can actually fill channel_models from other.channel_models.
    for (unsigned int c = 0; c < other.channel_models.size(); c++) {
        channel_models.push_back(new ChannelModel(*other.channel_models[c]));
    }
    return *this;
}

SequencingModel::SequencingModel(SequencingModel&& other) {
    p_edman_failure = move(other.p_edman_failure);
    p_detach = move(other.p_detach);
    p_initial_break_n = move(other.p_initial_break_n);
    p_cyclic_break_n = move(other.p_cyclic_break_n);
    channel_models = move(other.channel_models);
}

SequencingModel::~SequencingModel() {
    for (ChannelModel* channel_model : channel_models) {
        if (channel_model != NULL) {
            delete channel_model;
        }
    }
}

SequencingModel SequencingModel::with_mu_as_one() const {
    SequencingModel x;
    x.p_detach = p_detach;
    x.p_edman_failure = p_edman_failure;
    for (unsigned int i = 0; i < channel_models.size(); i++) {
        x.channel_models.push_back(
                new ChannelModel(channel_models[i]->with_mu_as_one()));
    }
    return x;
}

double SequencingModel::distance(
        const SequencingModel& sequencing_model) const {
    double dist = 0.0;
    dist = max(dist, abs(p_edman_failure - sequencing_model.p_edman_failure));
    dist = max(dist, abs(p_detach - sequencing_model.p_detach));
    dist = max(dist,
               abs(p_initial_break_n - sequencing_model.p_initial_break_n));
    dist = max(dist, abs(p_cyclic_break_n - sequencing_model.p_cyclic_break_n));
    for (unsigned int i = 0; i < channel_models.size(); i++) {
        dist = max(dist,
                   channel_models[i]->distance(
                           *sequencing_model.channel_models[i]));
    }
    return dist;
}

string SequencingModel::debug_string() const {
    string s = "";
    s += "Edman failure rate: " + to_string(p_edman_failure) + ", ";
    s += "Detach rate: " + to_string(p_detach) + ", ";
    s += "Initial blocked N-terminus rate: " + to_string(p_initial_break_n)
         + ", ";
    s += "Cyclic blocked N-terminus rate: " + to_string(p_cyclic_break_n);
    for (unsigned int i = 0; i < channel_models.size(); i++) {
        s += ", ";
        s += "Channel " + to_string(i) + " info: (";
        s += channel_models[i]->debug_string() + ")";
    }
    return s;
}

}  // namespace whatprot
