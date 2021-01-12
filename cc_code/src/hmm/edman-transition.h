/******************************************************************************\
* Author: Matthew Beauregard Smith                                             *
* Affiliation: The University of Texas at Austin                               *
* Department: Oden Institute and Institute for Cellular and Molecular Biology  *
* PI: Edward Marcotte                                                          *
* Project: Protein Fluorosequencing                                            *
\******************************************************************************/

#ifndef FLUOROSEQ_FWD_ALG_EDMAN_TRANSITION_H
#define FLUOROSEQ_FWD_ALG_EDMAN_TRANSITION_H

// Local project headers:
#include "common/dye-track.h"
#include "tensor/tensor.h"

namespace fluoroseq {

class EdmanTransition {
public:
    EdmanTransition(double p_edman_failure,
                    const DyeSeq& dye_seq,
                    const DyeTrack& dye_track);
    void forward(Tensor* tensor, int timestep) const;

    DyeSeq dye_seq;
    DyeTrack dye_track;
    double p_edman_failure;
};

}  // namespace fluoroseq

#endif  // FLUOROSEQ_FWD_ALG_EDMAN_TRANSITION_H