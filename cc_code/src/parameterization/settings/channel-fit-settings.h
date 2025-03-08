/******************************************************************************\
* Author: Matthew Beauregard Smith                                             *
* Affiliation: Erisyon, Inc.                                                   *
* Project: Protein Fluorosequencing                                            *
\******************************************************************************/

#ifndef WHATPROT_PARAMETERIZATION_SETTINGS_CHANNEL_FIT_SETTINGS_H
#define WHATPROT_PARAMETERIZATION_SETTINGS_CHANNEL_FIT_SETTINGS_H

namespace whatprot {

class ChannelFitSettings {
public:
    ChannelFitSettings();

    bool hold_p_bleach;
    bool hold_p_dud;
    // We assume that users will in all cases hold constant bg_sig, mu, and sig.
    // This could perhaps be revisited if the model is made to better match the
    // data or the data is made more clean.
};

}  // namespace whatprot

#endif  // WHATPROT_PARAMETERIZATION_SETTINGS_CHANNEL_FIT_SETTINGS_H