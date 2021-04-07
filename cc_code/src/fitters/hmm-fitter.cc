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
#include <iostream>
#include <vector>

// Local project headers:
#include "common/dye-seq.h"
#include "common/dye-track.h"
#include "common/error-model.h"
#include "common/radiometry.h"
#include "hmm/fit/error-model-fitter.h"
#include "hmm/hmm/peptide-hmm.h"
#include "hmm/hmm/stuck-dye-hmm.h"
#include "hmm/precomputations/dye-seq-precomputations.h"
#include "hmm/precomputations/radiometry-precomputations.h"
#include "hmm/precomputations/universal-precomputations.h"

namespace whatprot {

namespace {
using std::cout;
using std::move;
using std::vector;
}  // namespace

HMMFitter::HMMFitter(int num_timesteps,
                     int num_channels,
                     double stopping_threshold,
                     const ErrorModel& error_model,
                     const DyeSeq& dye_seq)
        : num_timesteps(num_timesteps),
          num_channels(num_channels),
          stopping_threshold(stopping_threshold),
          error_model(error_model),
          dye_seq(dye_seq) {
    max_num_dyes = 0;
    for (int c = 0; c < num_channels; c++) {
        int num_dyes = 0;
        for (int i = 0; i < dye_seq.length; i++) {
            if (dye_seq[i] == c) {
                num_dyes++;
            }
        }
        if (num_dyes > max_num_dyes) {
            max_num_dyes = num_dyes;
        }
    }
}

ErrorModel HMMFitter::fit(const std::vector<Radiometry>& radiometries) const {
    ErrorModel em = error_model;
    cout << em.debug_string() << "\n";
    while (true) {
        ErrorModelFitter fitter;
        DyeSeqPrecomputations dye_seq_precomputations(
                dye_seq, em, num_timesteps, num_channels);
        UniversalPrecomputations universal_precomputations(em, num_channels);
        universal_precomputations.set_max_num_dyes(max_num_dyes);
        for (const Radiometry& radiometry : radiometries) {
            RadiometryPrecomputations radiometry_precomputations(
                    radiometry, em, max_num_dyes);
            PeptideHMM hmm(num_timesteps,
                           num_channels,
                           dye_seq_precomputations,
                           radiometry_precomputations,
                           universal_precomputations);
            ErrorModelFitter peptide_fitter;
            double peptide_prob =
                    hmm.improve_fit(&peptide_fitter) * (1 - em.stuck_dye_ratio);
            vector<ErrorModelFitter> stuck_dye_fitters;
            vector<double> stuck_dye_probs;
            double total_prob = peptide_prob;
            for (int c = 0; c < num_channels; c++) {
                stuck_dye_fitters.push_back(ErrorModelFitter());
                StuckDyeHMM stuck_dye_hmm(num_timesteps,
                                          num_channels,
                                          c,
                                          radiometry_precomputations,
                                          universal_precomputations);
                stuck_dye_probs.push_back(
                        stuck_dye_hmm.improve_fit(&stuck_dye_fitters.back())
                        * em.stuck_dye_ratio / (double)num_channels);
                total_prob += stuck_dye_probs.back();
            }
            if (total_prob != 0.0) {
                peptide_fitter *= (peptide_prob / total_prob);
                fitter += peptide_fitter;
                for (int c = 0; c < num_channels; c++) {
                    stuck_dye_fitters[c] *= (stuck_dye_probs[c] / total_prob);
                    fitter += stuck_dye_fitters[c];
                }
                fitter.stuck_dye_ratio_fit.numerator += (1.0 - peptide_prob / total_prob);
                fitter.stuck_dye_ratio_fit.denominator += 1.0;
            }
        }
        ErrorModel next = fitter.error_model();
        cout << next.debug_string() << "\n";
        double relative_distance = em.relative_distance(next);
        cout << "relative distance: " << relative_distance << "\n";
        if (relative_distance < stopping_threshold) {
            return next;
        }
        em = next;
    }
}

}  // namespace whatprot
