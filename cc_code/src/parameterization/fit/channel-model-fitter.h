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

namespace whatprot {

class ChannelModelFitter {
public:
    ChannelModelFitter();
    ChannelModelFitter(const ChannelModelFitter& other);
    ChannelModelFitter(ChannelModelFitter&& other);
    ~ChannelModelFitter();
    ChannelModel get() const;
    ChannelModelFitter operator+(const ChannelModelFitter& other) const;
    void operator+=(const ChannelModelFitter& other);
    void operator*=(double weight_adjustment);
    ParameterFitter p_initial_bleach_fit;
    ParameterFitter p_cyclic_bleach_fit;
    ParameterFitter p_dud_fit;
    NormalDistributionFitter* distribution_fit;
};

}  // namespace whatprot

#endif  // WHATPROT_PARAMETERIZATION_FIT_CHANNEL_MODEL_FITTER_H
