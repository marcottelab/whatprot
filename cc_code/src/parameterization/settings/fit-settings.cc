/******************************************************************************\
* Author: Matthew Beauregard Smith                                             *
* Affiliation: Erisyon, Inc.                                                   *
* Project: Protein Fluorosequencing                                            *
\******************************************************************************/

// Defining symbols from header:
#include "fit-settings.h"

// Standard C++ library headers:
#include <fstream>
#include <string>
#include <utility>

// External headers:
#include "json.hpp"

// Local project headers:
#include "parameterization/settings/channel-fit-settings.h"

namespace {
using json = nlohmann::json;
using std::ifstream;
using std::move;
using std::string;
}  // namespace

namespace whatprot {

FitSettings::FitSettings() {}

FitSettings::FitSettings(unsigned int num_channels)
        : hold_p_edman_failure(false),
          hold_p_detach(false),
          hold_p_initial_detach(false),
          hold_p_initial_detach_decay(false),
          hold_p_initial_block(false),
          hold_p_cyclic_block(false) {
    for (unsigned int c = 0; c < num_channels; c++) {
        channel_fit_settings.push_back(new ChannelFitSettings());
    }
}

FitSettings::FitSettings(unsigned int num_channels,
                         const string& fit_settings_filename) {
    ifstream f(fit_settings_filename);
    json data = json::parse(f);
    if (data.contains("hold_p_edman_failure")) {
        hold_p_edman_failure = data["hold_p_edman_failure"].get<bool>();
    } else {
        hold_p_edman_failure = false;
    }
    if (data.contains("hold_p_detach")) {
        hold_p_detach = data["hold_p_detach"].get<bool>();
    } else {
        hold_p_detach = false;
    }
    if (data.contains("hold_p_initial_detach")) {
        hold_p_initial_detach = data["hold_p_initial_detach"].get<bool>();
    } else {
        hold_p_initial_detach = false;
    }
    if (data.contains("hold_p_initial_detach_decay")) {
        hold_p_initial_detach_decay =
                data["hold_p_initial_detach_decay"].get<bool>();
    } else {
        hold_p_initial_detach_decay = false;
    }
    if (data.contains("hold_p_initial_block")) {
        hold_p_initial_block = data["hold_p_initial_block"].get<bool>();
    } else {
        hold_p_initial_block = false;
    }
    if (data.contains("hold_p_cyclic_block")) {
        hold_p_cyclic_block = data["hold_p_cyclic_block"].get<bool>();
    } else {
        hold_p_cyclic_block = false;
    }
    if (data.contains("channel_settings")) {
        for (auto& channel_data : data["channel_settings"]) {
            channel_fit_settings.push_back(new ChannelFitSettings());
            ChannelFitSettings* cfsb = channel_fit_settings.back();
            if (channel_data.contains("hold_p_bleach")) {
                cfsb->hold_p_bleach = channel_data["hold_p_bleach"].get<bool>();
            } else {
                cfsb->hold_p_bleach = false;
            }
            if (channel_data.contains("hold_p_dud")) {
                cfsb->hold_p_dud = channel_data["hold_p_dud"].get<bool>();
            } else {
                cfsb->hold_p_dud = false;
            }
            // bg_sig, mu, and sig are never updated and don't have an
            // associated setting. See comment in channel-fit-settings.h.
        }
    } else {
        for (unsigned int c = 0; c < num_channels; c++) {
            channel_fit_settings.push_back(new ChannelFitSettings());
        }
    }
}

FitSettings::FitSettings(const FitSettings& other) {
    hold_p_edman_failure = other.hold_p_edman_failure;
    hold_p_detach = other.hold_p_detach;
    hold_p_initial_detach = other.hold_p_initial_detach;
    hold_p_initial_detach_decay = other.hold_p_initial_detach_decay;
    hold_p_initial_block = other.hold_p_initial_block;
    hold_p_cyclic_block = other.hold_p_cyclic_block;
    for (unsigned int c = 0; c < other.channel_fit_settings.size(); c++) {
        channel_fit_settings.push_back(
                new ChannelFitSettings(*other.channel_fit_settings[c]));
    }
}

FitSettings& FitSettings::operator=(const FitSettings& other) {
    hold_p_edman_failure = other.hold_p_edman_failure;
    hold_p_detach = other.hold_p_detach;
    hold_p_initial_detach = other.hold_p_initial_detach;
    hold_p_initial_detach_decay = other.hold_p_initial_detach_decay;
    hold_p_initial_block = other.hold_p_initial_block;
    hold_p_cyclic_block = other.hold_p_cyclic_block;
    // This function is not necessarily used as a constructor. It is very
    // important to clear contents of channel_fit_settings before filling it.
    for (ChannelFitSettings* cfs : channel_fit_settings) {
        if (cfs != NULL) {
            delete cfs;
        }
    }
    channel_fit_settings.resize(0);
    for (unsigned int c = 0; c < other.channel_fit_settings.size(); c++) {
        channel_fit_settings.push_back(
                new ChannelFitSettings(*other.channel_fit_settings[c]));
    }
    return *this;
}

FitSettings::FitSettings(FitSettings&& other) {
    hold_p_edman_failure = move(other.hold_p_edman_failure);
    hold_p_detach = move(other.hold_p_detach);
    hold_p_initial_detach = move(other.hold_p_initial_detach);
    hold_p_initial_detach_decay = move(other.hold_p_initial_detach_decay);
    hold_p_initial_block = move(other.hold_p_initial_block);
    hold_p_cyclic_block = move(other.hold_p_cyclic_block);
    channel_fit_settings = move(other.channel_fit_settings);
}

FitSettings::~FitSettings() {
    for (ChannelFitSettings* c_fs : channel_fit_settings) {
        if (c_fs != NULL) {
            delete c_fs;
        }
    }
}

}  // namespace whatprot
