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
#include "hmm/fit/log-normal-distribution-fitter.h"
#include "hmm/fit/normal-distribution-fitter.h"

namespace whatprot {

ErrorModelFitter::ErrorModelFitter() {
    distribution_fit = new LogNormalDistributionFitter();
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
                      distribution_fit->get_sigma(),
                      stuck_dye_ratio_fit.get(),
                      p_stuck_dye_loss_fit.get());
}

ErrorModelFitter ErrorModelFitter::operator+(
        const ErrorModelFitter& other) const {
    ErrorModelFitter result_fitter;
    result_fitter.p_edman_failure_fit =
            p_edman_failure_fit + other.p_edman_failure_fit;
    result_fitter.p_detach_fit = p_detach_fit + other.p_detach_fit;
    result_fitter.p_bleach_fit = p_bleach_fit + other.p_bleach_fit;
    result_fitter.p_dud_fit = p_dud_fit + other.p_dud_fit;
    *result_fitter.distribution_fit =
            *distribution_fit + *other.distribution_fit;
    result_fitter.stuck_dye_ratio_fit =
            stuck_dye_ratio_fit + other.stuck_dye_ratio_fit;
    result_fitter.p_stuck_dye_loss_fit =
            p_stuck_dye_loss_fit + other.p_stuck_dye_loss_fit;
    return result_fitter;
}

void ErrorModelFitter::operator+=(const ErrorModelFitter& other) {
    p_edman_failure_fit += other.p_edman_failure_fit;
    p_detach_fit += other.p_detach_fit;
    p_bleach_fit += other.p_bleach_fit;
    p_dud_fit += other.p_dud_fit;
    *distribution_fit += *other.distribution_fit;
    stuck_dye_ratio_fit += other.stuck_dye_ratio_fit;
    p_stuck_dye_loss_fit += other.p_stuck_dye_loss_fit;
}

void ErrorModelFitter::operator*=(double weight_adjustment) {
    p_edman_failure_fit *= weight_adjustment;
    p_detach_fit *= weight_adjustment;
    p_bleach_fit *= weight_adjustment;
    p_dud_fit *= weight_adjustment;
    *distribution_fit *= weight_adjustment;
    stuck_dye_ratio_fit *= weight_adjustment;
    p_stuck_dye_loss_fit *= weight_adjustment;
}

}  // namespace whatprot
