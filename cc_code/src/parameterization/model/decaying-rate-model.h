/******************************************************************************\
* Author: Matthew Beauregard Smith                                             *
* Affiliation: Erisyon, Inc.                                                   *
* Project: Protein Fluorosequencing                                            *
\******************************************************************************/

#ifndef WHATPROT_PARAMETERIZATION_MODEL_DECAYING_RATE_MODEL_H
#define WHATPROT_PARAMETERIZATION_MODEL_DECAYING_RATE_MODEL_H

// Standard C++ library headers:
#include <string>
#include <vector>

// Local project headers:
#include "parameterization/model/channel-model.h"

namespace whatprot {

class DecayingRateModel {
public:
    DecayingRateModel();
    DecayingRateModel(const DecayingRateModel& other);
    DecayingRateModel& operator=(const DecayingRateModel& other);
    DecayingRateModel(DecayingRateModel&& other);
    double operator[](unsigned int i) const;
    double distance(const DecayingRateModel& other) const;

    double base;
    double initial;
    double initial_decay;
};

}  // namespace whatprot

#endif  // WHATPROT_PARAMETERIZATION_MODEL_DECAYING_RATE_MODEL_H
