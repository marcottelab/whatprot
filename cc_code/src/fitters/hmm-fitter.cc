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
#include "parameterization/settings/channel-fit-settings.h"
#include "parameterization/settings/fit-settings.h"
#include "parameterization/settings/sequencing-settings.h"

namespace whatprot {

HMMFitter::HMMFitter(unsigned int num_timesteps,
                     unsigned int num_channels,
                     double stopping_threshold,
                     double max_runtime,
                     const SequencingModel& seq_model,
                     const SequencingSettings& seq_settings,
                     const FitSettings& fit_settings,
                     const DyeSeq& dye_seq)
        : dye_seq(dye_seq),
          seq_model(seq_model),
          seq_settings(seq_settings),
          fit_settings(fit_settings),
          stopping_threshold(stopping_threshold),
          max_runtime(max_runtime),
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

void HMMFitter::update_with_holds(const SequencingModel& update,
                                  SequencingModel* sm) const {
    if (!fit_settings.hold_p_edman_failure) {
        sm->p_edman_failure = update.p_edman_failure;
    }
    if (!fit_settings.hold_p_detach) {
        sm->p_detach = update.p_detach;
    }
    if (!fit_settings.hold_p_initial_block) {
        sm->p_initial_block = update.p_initial_block;
    }
    if (!fit_settings.hold_p_cyclic_block) {
        sm->p_cyclic_block = update.p_cyclic_block;
    }
    for (unsigned int i = 0; i < fit_settings.channel_fit_settings.size();
         i++) {
        const ChannelFitSettings& c_fs = *fit_settings.channel_fit_settings[i];
        const ChannelModel& c_update = *update.channel_models[i];
        ChannelModel* c_sm = sm->channel_models[i];
        if (!c_fs.hold_p_bleach) {
            c_sm->p_bleach = c_update.p_bleach;
        }
        if (!c_fs.hold_p_dud) {
            c_sm->p_dud = c_update.p_dud;
        }
        // bg_sig, mu, and sig are never updated. See comment in
        // channel-fit-settings.h.
    }
}

}  // namespace whatprot
