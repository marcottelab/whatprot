/******************************************************************************\
* Author: Matthew Beauregard Smith                                             *
* Affiliation: The University of Texas at Austin                               *
* Department: Oden Institute and Institute for Cellular and Molecular Biology  *
* PI: Edward Marcotte                                                          *
* Project: Protein Fluorosequencing                                            *
\******************************************************************************/

#ifndef WHATPROT_FWD_ALG_EDMAN_TRANSITION_H
#define WHATPROT_FWD_ALG_EDMAN_TRANSITION_H

// Local project headers:
#include "common/dye-track.h"
#include "tensor/tensor.h"

namespace whatprot {

class EdmanTransition {
public:
    EdmanTransition(double p_edman_failure,
                    const DyeSeq& dye_seq,
                    const DyeTrack& dye_track);
    void operator()(Tensor* tensor, int timestep) const;

    DyeSeq dye_seq;
    DyeTrack dye_track;
    double p_edman_failure;
};

}  // namespace whatprot

#endif  // WHATPROT_FWD_ALG_EDMAN_TRANSITION_H