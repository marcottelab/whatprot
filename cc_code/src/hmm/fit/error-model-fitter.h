/******************************************************************************\
* Author: Matthew Beauregard Smith                                             *
* Affiliation: The University of Texas at Austin                               *
* Department: Oden Institute and Institute for Cellular and Molecular Biology  *
* PI: Edward Marcotte                                                          *
* Project: Protein Fluorosequencing                                            *
\******************************************************************************/

#ifndef WHATPROT_HMM_FIT_ERROR_MODEL_FITTER_H
#define WHATPROT_HMM_FIT_ERROR_MODEL_FITTER_H

// Local project headers:
#include "hmm/fit/log-normal-distribution-fitter.h"
#include "hmm/fit/parameter-fitter.h"

namespace whatprot {

class ErrorModelFitter {
public:
    ErrorModelFitter();
    ErrorModelFitter(const ErrorModelFitter& other);
    ErrorModelFitter(ErrorModelFitter&& other);
    ~ErrorModelFitter();
    ErrorModel error_model() const;
    ErrorModelFitter operator+(const ErrorModelFitter& other) const;
    void operator+=(const ErrorModelFitter& other);
    void operator*=(double weight_adjustment);
    ParameterFitter p_edman_failure_fit;
    ParameterFitter p_detach_fit;
    ParameterFitter p_bleach_fit;
    ParameterFitter p_dud_fit;
    LogNormalDistributionFitter* distribution_fit;
    ParameterFitter stuck_dye_ratio_fit;
    ParameterFitter p_stuck_dye_loss_fit;
};

}  // namespace whatprot

#endif  // WHATPROT_HMM_FIT_ERROR_MODEL_FITTER_H
