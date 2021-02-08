/******************************************************************************\
* Author: Matthew Beauregard Smith                                             *
* Affiliation: The University of Texas at Austin                               *
* Department: Oden Institute and Institute for Cellular and Molecular Biology  *
* PI: Edward Marcotte                                                          *
* Project: Protein Fluorosequencing                                            *
\******************************************************************************/

#ifndef WHATPROT_HMM_PRECOMPUTATIONS_DYE_SEQ_PRECOMPUTATIONS_H
#define WHATPROT_HMM_PRECOMPUTATIONS_DYE_SEQ_PRECOMPUTATIONS_H

// Standard C++ library headers:
#include <vector>

// Local project headers:
#include "common/dye-seq.h"
#include "common/error-model.h"
#include "hmm/step/edman-transition.h"

namespace whatprot {

class DyeSeqPrecomputations {
public:
    DyeSeqPrecomputations(const DyeSeq& dye_seq,
                          const ErrorModel& error_model,
                          int num_timesteps,
                          int num_channels);
    std::vector<int> tensor_shape;
    EdmanTransition edman_transition;
};

}  // namespace whatprot

#endif  // WHATPROT_HMM_PRECOMPUTATIONS_DYE_SEQ_PRECOMPUTATIONS_H
