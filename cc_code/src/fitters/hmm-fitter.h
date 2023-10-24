/******************************************************************************\
* Author: Matthew Beauregard Smith                                             *
* Affiliation: The University of Texas at Austin                               *
* Department: Oden Institute and Institute for Cellular and Molecular Biology  *
* PI: Edward Marcotte                                                          *
* Project: Protein Fluorosequencing                                            *
\******************************************************************************/

#ifndef WHATPROT_FITTERS_HMM_FITTER_H
#define WHATPROT_FITTERS_HMM_FITTER_H

// Standard C++ library headers:
#include <cmath>
#include <iomanip>
#include <iostream>
#include <vector>

// Local project headers:
#include "common/dye-seq.h"
#include "common/radiometry.h"
#include "hmm/hmm/peptide-hmm.h"
#include "hmm/precomputations/dye-seq-precomputations.h"
#include "hmm/precomputations/radiometry-precomputations.h"
#include "hmm/precomputations/universal-precomputations.h"
#include "parameterization/model/sequencing-model.h"
#include "parameterization/settings/fit-settings.h"
#include "parameterization/settings/sequencing-settings.h"
#include "util/schroedinger-pointers.h"
#include "util/time.h"

namespace whatprot {

class HMMFitter {
public:
    HMMFitter(unsigned int num_timesteps,
              unsigned int num_channels,
              double stopping_threshold,
              double max_runtime,
              const SequencingModel& seq_model,
              const SequencingSettings& seq_settings,
              const FitSettings& fit_settings,
              const DyeSeq& dye_seq);

    // helper function
    void update_with_holds(const SequencingModel& update,
                           SequencingModel* sm) const;

    // Note: R must be Radiometry type or Radiometry pointer type.
    template <class R>
    double fit(const std::vector<R>& radiometries,
               SequencingModel* x,
               double* step_size) const {
        SequencingModel sm = seq_model;
        double start_time = wall_time();
        while (true) {
            SequencingModelFitter fitter(num_channels);
            DyeSeqPrecomputations dye_seq_precomputations(
                    dye_seq, sm, num_timesteps, num_channels);
            UniversalPrecomputations universal_precomputations(sm,
                                                               num_channels);
            universal_precomputations.set_max_num_dyes(max_num_dyes);
            for (const auto& radiometry : radiometries) {
                RadiometryPrecomputations radiometry_precomputations(
                        dereference_if_pointer(radiometry),
                        sm,
                        seq_settings,
                        max_num_dyes);
                PeptideHMM hmm(num_timesteps,
                               num_channels,
                               dye_seq_precomputations,
                               radiometry_precomputations,
                               universal_precomputations);
                SequencingModelFitter peptide_fitter(num_channels);
                hmm.improve_fit(&peptide_fitter);
                fitter += peptide_fitter;
            }
            // Here we perform a correction to account for the peptides that
            // wouldn't be seen due to all fluorophores being duds. This fixes
            // bias in result for p_dud on all channels.
            // TODO: maybe make this cleaner - break into separate function?
            // TODO: also inefficient, maybe should be a precomputation somehow?
            double ratio_hidden = 1.0;
            for (unsigned int i = 0; i < dye_seq.length; i++) {
                if (dye_seq[i] != -1) {
                    ratio_hidden *= sm.channel_models[dye_seq[i]]->p_dud;
                }
            }
            double magic_ratio = 1.0 / (1.0 - ratio_hidden) - 1.0;
            double expected_hidden_count = magic_ratio * radiometries.size();
            // We have to account for expected hidden count for EACH fluorophore
            // so that they are additive (i.e., two fluorophores equals double
            // the effect on the fitter).
            for (unsigned int i = 0; i < dye_seq.length; i++) {
                if (dye_seq[i] != -1) {
                    fitter.channel_fits[dye_seq[i]]->p_dud_fit.numerator +=
                            expected_hidden_count;
                    fitter.channel_fits[dye_seq[i]]->p_dud_fit.denominator +=
                            expected_hidden_count;
                }
            }
            SequencingModel next = sm;
            update_with_holds(fitter.get(), &next);
            *step_size = sm.distance(next);
            if (*step_size < stopping_threshold) {
                *x = next;
                break;
            }
            if (wall_time() - start_time > max_runtime) {
                *x = next;
                break;
            }
            sm = next;
        }
        return log_likelihood(radiometries, *x);
    }

    // Note: R must be Radiometry type or Radiometry pointer type.
    template <class R>
    double log_likelihood(const std::vector<R>& radiometries,
                          const SequencingModel& seq_model) const {
        DyeSeqPrecomputations dye_seq_precomputations(
                dye_seq, seq_model, num_timesteps, num_channels);
        UniversalPrecomputations universal_precomputations(seq_model,
                                                           num_channels);
        universal_precomputations.set_max_num_dyes(max_num_dyes);
        double log_l = 0.0;
        for (const auto& radiometry : radiometries) {
            RadiometryPrecomputations radiometry_precomputations(
                    dereference_if_pointer(radiometry),
                    seq_model,
                    seq_settings,
                    max_num_dyes);
            PeptideHMM hmm(num_timesteps,
                           num_channels,
                           dye_seq_precomputations,
                           radiometry_precomputations,
                           universal_precomputations);
            log_l += log(hmm.probability());
        }
        return log_l;
    }

    const DyeSeq& dye_seq;
    std::vector<DyeSeq> stuck_dyes;
    const SequencingModel& seq_model;
    const SequencingSettings& seq_settings;
    const FitSettings& fit_settings;
    double stopping_threshold;
    double max_runtime;
    unsigned int num_timesteps;
    unsigned int num_channels;
    int max_num_dyes;
};

}  // namespace whatprot

#endif  // WHATPROT_FITTERS_HMM_FITTER_H
