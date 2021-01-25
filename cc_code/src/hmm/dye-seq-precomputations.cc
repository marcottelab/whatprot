/******************************************************************************\
* Author: Matthew Beauregard Smith                                             *
* Affiliation: The University of Texas at Austin                               *
* Department: Oden Institute and Institute for Cellular and Molecular Biology  *
* PI: Edward Marcotte                                                          *
* Project: Protein Fluorosequencing                                            *
\******************************************************************************/

// Defining symbols from header:
#include "dye-seq-precomputations.h"

// Local project headers:
#include "common/dye-seq.h"
#include "common/dye-track.h"
#include "common/error-model.h"
#include "hmm/edman-transition.h"

namespace fluoroseq {

DyeSeqPrecomputations::DyeSeqPrecomputations(const DyeSeq& dye_seq,
                                             const ErrorModel& error_model,
                                             int num_timesteps,
                                             int num_channels)
        : edman_transition(error_model.p_edman_failure,
                           dye_seq,
                           DyeTrack(num_timesteps, num_channels, dye_seq)) {
    tensor_shape.resize(1 + num_channels);
    tensor_shape[0] = num_timesteps;
    for (int c = 0; c < num_channels; c++) {
        // This next line of code is a little confusing.
        //   * The first dimension of the tensor shape is always the timestep,
        //     so we need to add one to the channel to index to the correct
        //     dimension.
        //   * The zeroth step for the dye track gives us the maximum number of
        //     dyes possible, but the tensor shape for that channel needs to be
        //     one bigger than that to handle all values inclusively from 0 to
        //     the number of dyes.
        tensor_shape[1 + c] = 1 + edman_transition.dye_track(0, c);
    }
}

}  // namespace fluoroseq
