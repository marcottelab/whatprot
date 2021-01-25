/******************************************************************************\
* Author: Matthew Beauregard Smith                                             *
* Affiliation: The University of Texas at Austin                               *
* Department: Oden Institute and Institute for Cellular and Molecular Biology  *
* PI: Edward Marcotte                                                          *
* Project: Protein Fluorosequencing                                            *
\******************************************************************************/

#ifndef FLUOROSEQ_HMM_DYE_SEQ_PRECOMPUTATIONS_H
#define FLUOROSEQ_HMM_DYE_SEQ_PRECOMPUTATIONS_H

// Standard C++ library headers:
#include <vector>

// Local project headers:
#include "common/dye-seq.h"
#include "common/error-model.h"
#include "hmm/edman-transition.h"

namespace fluoroseq {

class DyeSeqPrecomputations {
public:
    DyeSeqPrecomputations(const DyeSeq& dye_seq,
                          const ErrorModel& error_model,
                          int num_timesteps,
                          int num_channels);
    std::vector<int> tensor_shape;
    EdmanTransition edman_transition;
};

}  // namespace fluoroseq

#endif  // FLUOROSEQ_HMM_DYE_SEQ_PRECOMPUTATIONS_H
