/******************************************************************************\
* Author: Matthew Beauregard Smith                                             *
* Affiliation: Erisyon, Inc.                                                   *
* Project: Protein Fluorosequencing                                            *
\******************************************************************************/

#ifndef WHATPROT_PARAMETERIZATION_SETTINGS_FIT_SETTINGS_H
#define WHATPROT_PARAMETERIZATION_SETTINGS_FIT_SETTINGS_H

// Standard C++ library headers:
#include <string>
#include <vector>

// Local project headers:
#include "parameterization/settings/channel-fit-settings.h"

namespace whatprot {

class FitSettings {
public:
    FitSettings();
    FitSettings(unsigned int num_channels);
    FitSettings(unsigned int num_channels,
                const std::string& fit_settings_filename);
    FitSettings(const FitSettings& other);
    FitSettings& operator=(const FitSettings& other);
    FitSettings(FitSettings&& other);
    ~FitSettings();

    bool hold_p_edman_failure;
    bool hold_p_detach;
    bool hold_p_initial_detach;
    bool hold_p_initial_detach_decay;
    bool hold_p_initial_block;
    bool hold_p_cyclic_block;
    std::vector<ChannelFitSettings*> channel_fit_settings;
};

}  // namespace whatprot

#endif  // WHATPROT_PARAMETERIZATION_SETTINGS_FIT_SETTINGS_H