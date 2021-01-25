/******************************************************************************\
* Author: Matthew Beauregard Smith                                             *
* Affiliation: The University of Texas at Austin                               *
* Department: Oden Institute and Institute for Cellular and Molecular Biology  *
* PI: Edward Marcotte                                                          *
* Project: Protein Fluorosequencing                                            *
\******************************************************************************/

#ifndef FLUOROSEQ_HMM_HMM_H
#define FLUOROSEQ_HMM_HMM_H

// Standard C++ library headers:
#include <vector>

// Local project headers:
#include "hmm/dye-seq-precomputations.h"
#include "hmm/radiometry-precomputations.h"
#include "hmm/step.h"
#include "hmm/universal-precomputations.h"

namespace fluoroseq {

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
    double probability();
    std::vector<int> tensor_shape;
    std::vector<const Step*> steps;
};

}  // namespace fluoroseq

#endif  // FLUOROSEQ_HMM_HMM_H
