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
#include "hmm/state-vector/stuck-dye-state-vector.h"
#include "parameterization/fit/sequencing-model-fitter.h"

namespace whatprot {

StuckDyeTransition::StuckDyeTransition(double loss_rate, int channel)
        : loss_rate(loss_rate), channel(channel) {}

void StuckDyeTransition::forward(int* num_edmans,
                                 StuckDyeStateVector* sdsv) const {
    sdsv->no_dye += sdsv->dye * loss_rate;
    sdsv->dye *= (1 - loss_rate);
    (*num_edmans)++;
}

void StuckDyeTransition::backward(const StuckDyeStateVector& input,
                                  int* num_edmans,
                                  StuckDyeStateVector* output) const {
    (*num_edmans)--;
    output->dye = (1 - loss_rate) * input.dye + loss_rate * input.no_dye;
    output->no_dye = input.no_dye;
}

void StuckDyeTransition::improve_fit(
        const StuckDyeStateVector& forward_sdsv,
        const StuckDyeStateVector& backward_sdsv,
        const StuckDyeStateVector& next_backward_sdsv,
        int num_edmans,
        double probability,
        SequencingModelFitter* fitter) const {
    fitter->channel_fits[channel]->p_stuck_dye_loss_fit.numerator +=
            forward_sdsv.dye * loss_rate * next_backward_sdsv.no_dye
            / probability;
    fitter->channel_fits[channel]->p_stuck_dye_loss_fit.denominator +=
            forward_sdsv.dye * backward_sdsv.dye / probability;
}

}  // namespace whatprot
