/******************************************************************************\
* Author: Matthew Beauregard Smith                                             *
* Affiliation: The University of Texas at Austin                               *
* Department: Oden Institute and Institute for Cellular and Molecular Biology  *
* PI: Edward Marcotte                                                          *
* Project: Protein Fluorosequencing                                            *
\******************************************************************************/

// Defining symbols from header:
#include "hmm-fitter.h"

// Standard C++ library headers:

// Local project headers:
#include "common/dye-seq.h"
#include "common/dye-track.h"
#include "common/radiometry.h"
#include "hmm/hmm/peptide-hmm.h"
#include "hmm/precomputations/dye-seq-precomputations.h"
#include "hmm/precomputations/radiometry-precomputations.h"
#include "hmm/precomputations/universal-precomputations.h"
#include "parameterization/fit/sequencing-model-fitter.h"
#include "parameterization/model/sequencing-model.h"
#include "parameterization/settings/sequencing-settings.h"

namespace whatprot {

HMMFitter::HMMFitter(unsigned int num_timesteps,
                     unsigned int num_channels,
                     double stopping_threshold,
                     const SequencingModel& seq_model,
                     const SequencingSettings& seq_settings,
                     const DyeSeq& dye_seq)
        : dye_seq(dye_seq),
          seq_model(seq_model),
          seq_settings(seq_settings),
          stopping_threshold(stopping_threshold),
          num_timesteps(num_timesteps),
          num_channels(num_channels) {
    max_num_dyes = 0;
    for (unsigned int c = 0; c < num_channels; c++) {
        int num_dyes = 0;
        for (unsigned int i = 0; i < dye_seq.length; i++) {
            // Can safely compare with type-cast. c is unsigned because it is
            // a channel index, the DyeSeq [] operator is signed because values
            // of -1 indicate no dye in that position.
            if (dye_seq[i] == (int)c) {
                num_dyes++;
            }
        }
        if (num_dyes > max_num_dyes) {
            max_num_dyes = num_dyes;
        }
    }
}

}  // namespace whatprot
