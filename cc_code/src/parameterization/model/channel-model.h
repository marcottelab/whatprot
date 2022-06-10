/******************************************************************************\
* Author: Matthew Beauregard Smith                                             *
* Affiliation: The University of Texas at Austin                               *
* Department: Oden Institute and Institute for Cellular and Molecular Biology  *
* PI: Edward Marcotte                                                          *
* Project: Protein Fluorosequencing                                            *
\******************************************************************************/

#ifndef WHATPROT_PARAMETERIZATION_MODEL_CHANNEL_MODEL_H
#define WHATPROT_PARAMETERIZATION_MODEL_CHANNEL_MODEL_H

// Standard C++ library headers:
#include <functional>
#include <string>

namespace whatprot {

class ChannelModel {
public:
    virtual ~ChannelModel();
    virtual double pdf(double observed, int state) const;
    virtual double sigma(int state) const;
    double relative_distance(const ChannelModel& channel_model) const;
    std::string debug_string() const;

    double p_bleach;
    double p_dud;
    double bg_sig;
    double mu;
    double sig;
};

}  // namespace whatprot

#endif  // WHATPROT_PARAMETERIZATION_MODEL_CHANNEL_MODEL_H
