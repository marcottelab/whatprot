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
#include "hmm/step/edman-transition.h"
#include "parameterization/model/sequencing-model.h"

namespace whatprot {

class DyeSeqPrecomputations {
public:
    DyeSeqPrecomputations(const DyeSeq& dye_seq,
                          const SequencingModel& seq_model,
                          unsigned int num_timesteps,
                          unsigned int num_channels);
    // We forbid copy construction to guarantee consistent location of dye_track
    // for the edman_transition we are creating.
    DyeSeqPrecomputations(const DyeSeqPrecomputations& other) = delete;
    std::vector<unsigned int> tensor_shape;
    DyeTrack dye_track;  // MUST be before edman_transition
    EdmanTransition edman_transition;
};

}  // namespace whatprot

#endif  // WHATPROT_HMM_PRECOMPUTATIONS_DYE_SEQ_PRECOMPUTATIONS_H
