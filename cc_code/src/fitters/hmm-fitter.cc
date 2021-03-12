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
#include "hmm/hmm.h"
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
    // This is ugly. Modeling dyes stuck to the slide by creating a peptide
    // with more undyed amino acids than there are Edman Cycles. May need to
    // rework this for higher quality results.
    for (int c = 0; c < num_channels; c++) {
        char* dyes = new char[num_timesteps];
        for (int i = 0; i < num_timesteps - 1; i++) {
            dyes[i] = '.';
        }
        // This is a dirty trick to quickly turn c from an int to its
        // corresponding char.
        dyes[num_timesteps - 1] = '0' + c;
        stuck_dyes.emplace_back(num_timesteps, dyes);
        delete[] dyes;
    }
}

ErrorModel HMMFitter::fit(const std::vector<Radiometry>& radiometries) const {
    ErrorModel em = error_model;
    ErrorModel em_stuck_dye = error_model;
    em_stuck_dye.stuck_dye_ratio = 1.0;
    cout << em.debug_string() << "\n";
    while (true) {
        ErrorModelFitter fitter;
        ErrorModelFitter fitter_stuck_dye;
        DyeSeqPrecomputations dye_seq_precomputations(
                dye_seq, em, num_timesteps, num_channels);
        vector<DyeSeqPrecomputations> stuck_dye_precomputations;
        for (int c = 0; c < num_channels; c++) {
            stuck_dye_precomputations.emplace_back(
                    stuck_dyes[c], em_stuck_dye, num_timesteps, num_channels);
        }
        UniversalPrecomputations universal_precomputations(em, num_channels);
        UniversalPrecomputations universal_precomputations_stuck_dye(em_stuck_dye, num_channels);
        universal_precomputations.set_max_num_dyes(max_num_dyes);
        universal_precomputations_stuck_dye.set_max_num_dyes(max_num_dyes);
        for (const Radiometry& radiometry : radiometries) {
            ErrorModelFitter rad_fitter;
            ErrorModelFitter rad_fitter_stuck_dye;
            RadiometryPrecomputations radiometry_precomputations(
                    radiometry, em, max_num_dyes);
            RadiometryPrecomputations radiometry_precomputations_stuck_dye(
                    radiometry, em_stuck_dye, 1);
            HMM hmm(num_timesteps,
                    num_channels,
                    dye_seq_precomputations,
                    radiometry_precomputations,
                    universal_precomputations);
            double peptide_prob =
                    hmm.improve_fit(&rad_fitter) * (1 - em.stuck_dye_ratio);
            double total_prob = peptide_prob;
            for (int c = 0; c < num_channels; c++) {
                ErrorModelFitter temp_fitter;
                HMM stuck_dye_hmm(num_timesteps,
                                  num_channels,
                                  stuck_dye_precomputations[c],
                                  radiometry_precomputations_stuck_dye,
                                  universal_precomputations_stuck_dye);
                double stuck_dye_prob = stuck_dye_hmm.improve_fit(&temp_fitter)
                                         * em.stuck_dye_ratio
                                         / (double)num_channels;
                if (total_prob + stuck_dye_prob != 0.0) {
                    rad_fitter *= total_prob / (total_prob + stuck_dye_prob);
                    rad_fitter_stuck_dye *= total_prob / (total_prob + stuck_dye_prob);
                    temp_fitter *= stuck_dye_prob / (total_prob + stuck_dye_prob);
                }
                rad_fitter_stuck_dye += temp_fitter;
                total_prob += stuck_dye_prob;
            }
            if (total_prob != 0.0) {
                rad_fitter.stuck_dye_ratio_fit.numerator = (total_prob - peptide_prob) / total_prob;
                rad_fitter.stuck_dye_ratio_fit.denominator = 1.0;
            }
            fitter += rad_fitter;
            // fitter_stuck_dye += rad_fitter_stuck_dye;
            fitter += rad_fitter_stuck_dye;
        }
        ErrorModel next = fitter.error_model();
        // ErrorModel next_stuck_dye = fitter_stuck_dye.error_model();
        // next_stuck_dye.stuck_dye_ratio = 1.0;
        cout << next.debug_string() << "\n";
        // cout << next_stuck_dye.debug_string() << "\n";
        double relative_distance = em.relative_distance(next);
        // double relative_distance_stuck_dye = em_stuck_dye.relative_distance(next_stuck_dye);
        cout << "relative distance: " << relative_distance << "\n";
        // cout << "relative distance: " << relative_distance_stuck_dye << "\n";
        if (relative_distance < stopping_threshold ) { // && relative_distance_stuck_dye < stopping_threshold) {
            return next;
        }
        em = next;
        em_stuck_dye = em;
        // em_stuck_dye = next_stuck_dye;
    }
}

}  // namespace whatprot
