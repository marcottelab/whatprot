/******************************************************************************\
* Author: Matthew Beauregard Smith                                             *
* Affiliation: The University of Texas at Austin                               *
* Department: Oden Institute and Institute for Cellular and Molecular Biology  *
* PI: Edward Marcotte                                                          *
* Project: Protein Fluorosequencing                                            *
\******************************************************************************/

// Defining symbols from header:
#include "error-model-fitter.h"

// Local project headers:
#include "common/error-model.h"
#include "hmm/fit/distribution-fitter.h"
#include "hmm/fit/log-normal-distribution-fitter.h"
#include "hmm/fit/normal-distribution-fitter.h"

namespace fluoroseq {

ErrorModelFitter::ErrorModelFitter(DistributionType distribution_type) {
    switch (distribution_type) {
        case DistributionType::LOGNORMAL:
            distribution_fit = new LogNormalDistributionFitter();
            break;
        case DistributionType::NORMAL:
        default:
            distribution_fit = new NormalDistributionFitter();
            break;
    }
}

ErrorModelFitter::~ErrorModelFitter() {
    delete distribution_fit;
}

ErrorModel ErrorModelFitter::error_model() const {
    return ErrorModel(p_edman_failure_fit.get(),
                      p_detach_fit.get(),
                      p_bleach_fit.get(),
                      p_dud_fit.get(),
                      distribution_fit->get_type(),
                      distribution_fit->get_mu(),
                      distribution_fit->get_sigma());
}

}  // namespace fluoroseq
