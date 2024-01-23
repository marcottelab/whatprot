/******************************************************************************\
* Author: Matthew Beauregard Smith                                             *
* Affiliation: The University of Texas at Austin                               *
* Department: Oden Institute and Institute for Cellular and Molecular Biology  *
* PI: Edward Marcotte                                                          *
* Project: Protein Fluorosequencing                                            *
\******************************************************************************/

#ifndef WHATPROT_PARAMETERIZATION_FIT_CHANNEL_MODEL_FITTER_H
#define WHATPROT_PARAMETERIZATION_FIT_CHANNEL_MODEL_FITTER_H

// Local project headers:
#include "parameterization/fit/normal-distribution-fitter.h"
#include "parameterization/fit/parameter-fitter.h"
#include "parameterization/model/channel-model.h"

namespace whatprot {

class ChannelModelFitter {
public:
    ChannelModelFitter(const ChannelModel& prev);
    ChannelModelFitter(const ChannelModelFitter& other);
    ChannelModelFitter(ChannelModelFitter&& other);
    ChannelModel get() const;
    ChannelModelFitter operator+(const ChannelModelFitter& other) const;
    void operator+=(const ChannelModelFitter& other);
    void operator*=(double weight_adjustment);
    ParameterFitter p_bleach_fit;
    ParameterFitter p_dud_fit;
    ChannelModel prev;
};

}  // namespace whatprot

#endif  // WHATPROT_PARAMETERIZATION_FIT_CHANNEL_MODEL_FITTER_H
