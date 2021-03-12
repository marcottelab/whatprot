/******************************************************************************\
* Author: Matthew Beauregard Smith                                             *
* Affiliation: The University of Texas at Austin                               *
* Department: Oden Institute and Institute for Cellular and Molecular Biology  *
* PI: Edward Marcotte                                                          *
* Project: Protein Fluorosequencing                                            *
\******************************************************************************/

#ifndef WHATPROT_HMM_HMM_H
#define WHATPROT_HMM_HMM_H

// Standard C++ library headers:
#include <vector>

// Local project headers:
#include "hmm/fit/error-model-fitter.h"
#include "hmm/precomputations/dye-seq-precomputations.h"
#include "hmm/precomputations/radiometry-precomputations.h"
#include "hmm/precomputations/universal-precomputations.h"
#include "hmm/step/step.h"

namespace whatprot {

class HMM {
public:
    HMM(int num_timesteps,
        int num_channels,
        const DyeSeqPrecomputations& dye_seq_precomputations,
        const RadiometryPrecomputations& radiometry_precomputations,
        const UniversalPrecomputations& universal_precomputations);
    // This computes the probability of the provided dye seq producing the
    // provided radiometry. To do this efficiently, it uses a modified version
    // of the forward algorithm.
    double probability() const;
    // This will fit the data the HMM was provided with. It also computes the
    // probability as a side effect, so it returns this in case that is useful
    // to the caller.
    double improve_fit(ErrorModelFitter* fitter) const;
    std::vector<int> tensor_shape;
    std::vector<const Step*> steps;
};

}  // namespace whatprot

#endif  // WHATPROT_HMM_HMM_H
