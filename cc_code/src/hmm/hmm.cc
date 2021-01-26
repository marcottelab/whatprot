/******************************************************************************\
* Author: Matthew Beauregard Smith                                             *
* Affiliation: The University of Texas at Austin                               *
* Department: Oden Institute and Institute for Cellular and Molecular Biology  *
* PI: Edward Marcotte                                                          *
* Project: Protein Fluorosequencing                                            *
\******************************************************************************/

// Defining symbols from header:
#include "hmm.h"

// Standard C++ library headers:
#include <vector>

// Local project headers:
#include "hmm/step/binomial-transition.h"
#include "hmm/step/detach-transition.h"
#include "hmm/precomputations/dye-seq-precomputations.h"
#include "hmm/step/edman-transition.h"
#include "hmm/step/emission.h"
#include "hmm/precomputations/radiometry-precomputations.h"
#include "hmm/step/start.h"
#include "hmm/precomputations/universal-precomputations.h"
#include "tensor/tensor.h"

namespace fluoroseq {

namespace {
using std::vector;
}  // namespace

HMM::HMM(int num_timesteps,
         int num_channels,
         const DyeSeqPrecomputations& dye_seq_precomputations,
         const RadiometryPrecomputations& radiometry_precomputations,
         const UniversalPrecomputations& universal_precomputations) {
    steps.push_back(&universal_precomputations.start);
    for (int c = 0; c < num_channels; c++) {
        steps.push_back(&universal_precomputations.dud_transitions[c]);
    }
    steps.push_back(&radiometry_precomputations.emission);
    for (int t = 1; t < num_timesteps; t++) {
        steps.push_back(&universal_precomputations.detach_transition);
        for (int c = 0; c < num_channels; c++) {
            steps.push_back(&universal_precomputations.bleach_transitions[c]);
        }
        steps.push_back(&dye_seq_precomputations.edman_transition);
        steps.push_back(&radiometry_precomputations.emission);
    }
    tensor_shape = dye_seq_precomputations.tensor_shape;
}

double HMM::probability() {
    vector<const Step*>::iterator step = steps.begin();
    Tensor tensor(tensor_shape.size(), &tensor_shape[0]);
    int num_edmans = 0;
    while (step != steps.end()) {
        (*step)->forward(tensor, &num_edmans, &tensor);
        step++;
    }
    return tensor.sum();
}

}  // namespace fluoroseq
