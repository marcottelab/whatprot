/******************************************************************************\
* Author: Matthew Beauregard Smith                                             *
* Affiliation: The University of Texas at Austin                               *
* Department: Oden Institute and Institute for Cellular and Molecular Biology  *
* PI: Edward Marcotte                                                          *
* Project: Protein Fluorosequencing                                            *
\******************************************************************************/

#ifndef WHATPROT_HMM_HMM_PEPTIDE_HMM_H
#define WHATPROT_HMM_HMM_PEPTIDE_HMM_H

// Standard C++ library headers:
#include <vector>

// Local project headers:
#include "hmm/fit/error-model-fitter.h"
#include "hmm/hmm/generic-hmm.h"
#include "hmm/precomputations/dye-seq-precomputations.h"
#include "hmm/precomputations/radiometry-precomputations.h"
#include "hmm/precomputations/universal-precomputations.h"
#include "hmm/state-vector/peptide-state-vector.h"
#include "hmm/step/step.h"

namespace whatprot {

class PeptideHMM : public GenericHMM<PeptideStateVector> {
public:
    PeptideHMM(int num_timesteps,
               int num_channels,
               const DyeSeqPrecomputations& dye_seq_precomputations,
               const RadiometryPrecomputations& radiometry_precomputations,
               const UniversalPrecomputations& universal_precomputations);
    virtual PeptideStateVector create_states() const override;
    std::vector<int> tensor_shape;
};

}  // namespace whatprot

#endif  // WHATPROT_HMM_HMM_PEPTIDE_HMM_H
