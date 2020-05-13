// Author: Matthew Beauregard Smith
#ifndef FLUOROSEQ_FWD_ALG_EDMAN_TRANSITION_H
#define FLUOROSEQ_FWD_ALG_EDMAN_TRANSITION_H

#include "common/dye_track.h"
#include "tensor/tensor.h"

namespace fluoroseq {

class EdmanTransition {
public:
    EdmanTransition(double p_edman_failure,
                    const DyeSeq& dye_seq,
                    const DyeTrack& dye_track,
                    int max_failed_edmans);
    void operator()(Tensor* tensor, int timestep) const;

    DyeSeq dye_seq;
    DyeTrack dye_track;
    double p_edman_failure;
    int max_failed_edmans;
};

}  // namespace fluoroseq

#endif  // FLUOROSEQ_FWD_ALG_EDMAN_TRANSITION_H