/******************************************************************************\
* Author: Matthew Beauregard Smith                                             *
* Affiliation: The University of Texas at Austin                               *
* Department: Oden Institute and Institute for Cellular and Molecular Biology  *
* PI: Edward Marcotte                                                          *
* Project: Protein Fluorosequencing                                            *
\******************************************************************************/

#ifndef WHATPROT_PARAMETERIZATION_FIT_SEQUENCING_MODEL_FITTER_H
#define WHATPROT_PARAMETERIZATION_FIT_SEQUENCING_MODEL_FITTER_H

// Standard C++ library headers:
#include <vector>

// Local project headers:
#include "parameterization/fit/channel-model-fitter.h"
#include "parameterization/fit/log-normal-distribution-fitter.h"
#include "parameterization/fit/parameter-fitter.h"

namespace whatprot {

class SequencingModelFitter {
public:
    SequencingModelFitter();
    SequencingModelFitter(unsigned int num_channels);
    SequencingModelFitter(const SequencingModelFitter& other);
    SequencingModelFitter(SequencingModelFitter&& other);
    ~SequencingModelFitter();
    SequencingModel get() const;
    SequencingModelFitter operator+(const SequencingModelFitter& other) const;
    void operator+=(const SequencingModelFitter& other);
    void operator*=(double weight_adjustment);
    ParameterFitter p_edman_failure_fit;
    ParameterFitter p_initial_detach_fit;
    ParameterFitter p_cyclic_detach_fit;
    ParameterFitter p_initial_break_n_fit;
    ParameterFitter p_cyclic_break_n_fit;
    std::vector<ChannelModelFitter*> channel_fits;
};

}  // namespace whatprot

#endif  // WHATPROT_PARAMETERIZATION_FIT_SEQUENCING_MODEL_FITTER_H
