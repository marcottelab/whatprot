/******************************************************************************\
* Author: Matthew Beauregard Smith                                             *
* Affiliation: The University of Texas at Austin                               *
* Department: Oden Institute and Institute for Cellular and Molecular Biology  *
* PI: Edward Marcotte                                                          *
* Project: Protein Fluorosequencing                                            *
\******************************************************************************/

#ifndef FLUOROSEQ_HMM_FIT_ERROR_MODEL_FITTER_H
#define FLUOROSEQ_HMM_FIT_ERROR_MODEL_FITTER_H

// Local project headers:
#include "common/error-model.h"
#include "hmm/fit/distribution-fitter.h"
#include "hmm/fit/parameter-fitter.h"

namespace fluoroseq {

class ErrorModelFitter {
public:
    ErrorModelFitter(DistributionType distribution_type);
    ~ErrorModelFitter();
    ErrorModel error_model() const;
    ParameterFitter p_edman_failure_fit;
    ParameterFitter p_detach_fit;
    ParameterFitter p_bleach_fit;
    ParameterFitter p_dud_fit;
    DistributionFitter* distribution_fit;
};

}  // namespace fluoroseq

#endif  // FLUOROSEQ_HMM_FIT_ERROR_MODEL_FITTER_H
