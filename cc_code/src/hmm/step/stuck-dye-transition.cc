/******************************************************************************\
* Author: Matthew Beauregard Smith                                             *
* Affiliation: The University of Texas at Austin                               *
* Department: Oden Institute and Institute for Cellular and Molecular Biology  *
* PI: Edward Marcotte                                                          *
* Project: Protein Fluorosequencing                                            *
\******************************************************************************/

// Defining symbols from header:
#include "stuck-dye-transition.h"

// Local project headers:
#include "hmm/fit/error-model-fitter.h"
#include "hmm/state-vector/stuck-dye-state-vector.h"

namespace whatprot {

StuckDyeTransition::StuckDyeTransition(double loss_rate)
         : loss_rate(loss_rate) {}

void StuckDyeTransition::forward(StuckDyeStateVector* sdsv) const {
    sdsv->no_dye += sdsv->dye * loss_rate;
    sdsv->dye *= (1 - loss_rate);
    sdsv->num_edmans += 1;
}

void StuckDyeTransition::backward(const StuckDyeStateVector& input,
                                  StuckDyeStateVector* output) const {
    output->num_edmans = input.num_edmans - 1;
    output->dye = (1 - loss_rate) * input.dye + loss_rate * input.no_dye;
    output->no_dye = input.no_dye;
}

void StuckDyeTransition::improve_fit(const StuckDyeStateVector& forward_sdsv,
                                    const StuckDyeStateVector& backward_sdsv,
                                    const StuckDyeStateVector& next_backward_sdsv,
                                    double probability,
                                    ErrorModelFitter* fitter) const {
    fitter->p_stuck_dye_loss_fit.numerator += forward_sdsv.dye * loss_rate
                   * next_backward_sdsv.no_dye / probability;
    fitter->p_stuck_dye_loss_fit.denominator += forward_sdsv.dye * backward_sdsv.dye
                  / probability;
}

}  // namespace whatprot
