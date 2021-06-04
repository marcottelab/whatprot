/******************************************************************************\
* Author: Matthew Beauregard Smith                                             *
* Affiliation: The University of Texas at Austin                               *
* Department: Oden Institute and Institute for Cellular and Molecular Biology  *
* PI: Edward Marcotte                                                          *
* Project: Protein Fluorosequencing                                            *
\******************************************************************************/

#ifndef WHATPROT_HMM_STEP_STEP_H
#define WHATPROT_HMM_STEP_STEP_H

// Local project headers:
#include "parameterization/fit/sequencing-model-fitter.h"

namespace whatprot {

template <typename SV>  // SV is the state vector type.
class Step {
public:
    virtual void forward(int* num_edmans, SV* states) const = 0;
    virtual void backward(const SV& input,
                          int* num_edmans,
                          SV* output) const = 0;
    virtual void improve_fit(const SV& forward_states,
                             const SV& backward_states,
                             const SV& next_backward_states,
                             int num_edmans,
                             double probability,
                             SequencingModelFitter* fitter) const = 0;
};

}  // namespace whatprot

#endif  // WHATPROT_HMM_STEP_STEP_H
