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
        ErrorModelFitter fitter(DistributionType::LOGNORMAL);
        DyeSeqPrecomputations dye_seq_precomputations(
                dye_seq, em, num_timesteps, num_channels);
        UniversalPrecomputations universal_precomputations(em, num_channels);
        universal_precomputations.set_max_num_dyes(max_num_dyes);
        for (const Radiometry& radiometry : radiometries) {
            RadiometryPrecomputations radiometry_precomputations(
                    radiometry, em, max_num_dyes);
            HMM hmm(num_timesteps,
                    num_channels,
                    dye_seq_precomputations,
                    radiometry_precomputations,
                    universal_precomputations);
            hmm.improve_fit(&fitter);
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
