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
#include <vector>

// Local project headers:
#include "common/dye-seq.h"
#include "common/error-model.h"
#include "common/radiometry.h"
#include "hmm/hmm.h"
#include "hmm/precomputations/dye-seq-precomputations.h"
#include "hmm/precomputations/radiometry-precomputations.h"
#include "hmm/precomputations/universal-precomputations.h"

namespace whatprot {

class HMMFitter {
public:
    HMMFitter(int num_timesteps,
              int num_channels,
              double stopping_threshold,
              const ErrorModel& error_model,
              const DyeSeq& dye_seq);
    ErrorModel fit(const std::vector<Radiometry>& radiometries) const;
    const DyeSeq& dye_seq;
    const ErrorModel& error_model;
    double stopping_threshold;
    int num_timesteps;
    int num_channels;
    int max_num_dyes;
};

}  // namespace whatprot

#endif  // WHATPROT_FITTERS_HMM_FITTER_H
