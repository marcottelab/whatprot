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
#include <utility>

// External headers:
#include "json.hpp"

// Local project headers:
#include "parameterization/model/channel-model.h"
#include "parameterization/model/decaying-rate-model.h"

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

SequencingModel::SequencingModel(unsigned int num_channels) : p_detach() {
    for (unsigned int c = 0; c < num_channels; c++) {
        channel_models.push_back(new ChannelModel(c, num_channels));
        channel_models.back()->interactions.resize(num_channels, 1.0);
    }
}

SequencingModel::SequencingModel(const string& seq_model_filename) {
    ifstream f(seq_model_filename);
    json data = json::parse(f);
    p_edman_failure = data["p_edman_failure"].get<double>();
    p_detach.base = data["p_detach"].get<double>();
    p_detach.initial = data["p_initial_detach"].get<double>();
    p_detach.initial_decay = data["p_initial_detach_decay"].get<double>();
    p_initial_block = data["p_initial_block"].get<double>();
    p_cyclic_block = data["p_cyclic_block"].get<double>();
    unsigned int c = 0;
    unsigned int num_channels = data["channel_models"].size();
    // It would be cleaner to delegate this logic to ChannelModel, but it is not
    // clear what type "auto& channel_data" represents. The user-guide for
    // nlohmann::json provides no further insight. This works though, so it's
    // good enough.
    for (auto& channel_data : data["channel_models"]) {
        channel_models.push_back(new ChannelModel(c, num_channels));
        channel_models.back()->p_bleach =
                channel_data["p_bleach"].get<double>();
        channel_models.back()->p_dud = channel_data["p_dud"].get<double>();
        channel_models.back()->bg_sig = channel_data["bg_sig"].get<double>();
        channel_models.back()->mu = channel_data["mu"].get<double>();
        channel_models.back()->sig = channel_data["sig"].get<double>();
        channel_models.back()->interactions.resize(num_channels, 1.0);
        // Not required for users to fill in every possible interaction. The
        // default value for these is 0.0 (see channel-model.cc).
        for (auto& interaction : channel_data["interactions"]) {
            unsigned int other_channel = interaction["other_channel"];
            // Neither effect nor flat_effect need be specified. Optional. This
            // makes it easier to specify just one (would be weird to specify
            // neither but not a huge problem if someone does that).
            if (interaction.contains("effect")) {
                double effect = interaction["effect"];
                channel_models.back()->interactions[other_channel] =
                        1.0 - effect;
            }
            if (interaction.contains("flat_effect")) {
                double flat_effect = interaction["flat_effect"];
                channel_models.back()->flat_interactions[other_channel] =
                        1.0 - flat_effect;
            }
        }
        c++;
    }
}

SequencingModel::SequencingModel(const SequencingModel& other) {
    p_edman_failure = other.p_edman_failure;
    p_detach = other.p_detach;
    p_initial_block = other.p_initial_block;
    p_cyclic_block = other.p_cyclic_block;
    for (unsigned int c = 0; c < other.channel_models.size(); c++) {
        channel_models.push_back(new ChannelModel(*other.channel_models[c]));
    }
}

SequencingModel& SequencingModel::operator=(const SequencingModel& other) {
    p_edman_failure = other.p_edman_failure;
    p_detach = other.p_detach;
    p_initial_block = other.p_initial_block;
    p_cyclic_block = other.p_cyclic_block;
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
    p_initial_block = move(other.p_initial_block);
    p_cyclic_block = move(other.p_cyclic_block);
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
    x.p_initial_block = p_initial_block;
    x.p_cyclic_block = p_cyclic_block;
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
    dist = max(dist, p_detach.distance(sequencing_model.p_detach));
    dist = max(dist, abs(p_initial_block - sequencing_model.p_initial_block));
    dist = max(dist, abs(p_cyclic_block - sequencing_model.p_cyclic_block));
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
    s += "Detach rate: " + to_string(p_detach.base) + ", ";
    s += "Initial detach rate: " + to_string(p_detach.initial) + ", ";
    s += "Initial detach decay rate: " + to_string(p_detach.initial_decay)
         + ", ";
    s += "Initial blocked N-terminus rate: " + to_string(p_initial_block)
         + ", ";
    s += "Cyclic blocked N-terminus rate: " + to_string(p_cyclic_block);
    for (unsigned int i = 0; i < channel_models.size(); i++) {
        s += ", ";
        s += "Channel " + to_string(i) + " info: (";
        s += channel_models[i]->debug_string() + ")";
    }
    return s;
}

}  // namespace whatprot
