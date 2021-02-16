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
#include "hmm/fit/error-model-fitter.h"
#include "hmm/precomputations/dye-seq-precomputations.h"
#include "hmm/precomputations/radiometry-precomputations.h"
#include "hmm/precomputations/universal-precomputations.h"
#include "hmm/step/binomial-transition.h"
#include "hmm/step/detach-transition.h"
#include "hmm/step/edman-transition.h"
#include "hmm/step/emission.h"
#include "hmm/step/finish.h"
#include "hmm/step/start.h"
#include "tensor/tensor.h"

namespace whatprot {

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
    steps.push_back(&universal_precomputations.finish);
    tensor_shape = dye_seq_precomputations.tensor_shape;
}

double HMM::probability() const {
    vector<const Step*>::const_iterator step = steps.begin();
    Tensor tensor(tensor_shape.size(), &tensor_shape[0]);
    int num_edmans = 0;
    while (step != steps.end()) {
        (*step)->forward(tensor, &num_edmans, &tensor);
        step++;
    }
    return tensor.sum();
}

void HMM::improve_fit(ErrorModelFitter* fitter) const {
    vector<const Step*>::const_iterator step = steps.end();
    vector<Tensor> backward_tensors;
    backward_tensors.reserve(steps.size());
    // For efficiency, backwards_tensors is in the reverse order of what we
    // would like. Yes this is confusing...
    backward_tensors.emplace_back(tensor_shape.size(), &tensor_shape[0]);
    int num_edmans = tensor_shape[0] - 1;  // tensor_shape[0] is num_timesteps.
    while (step != steps.begin()) {
        step--;
        backward_tensors.emplace_back(tensor_shape.size(), &tensor_shape[0]);
        Tensor* left_tensor = &backward_tensors.back();
        const Tensor& right_tensor = *(&backward_tensors.back() - 1);
        (*step)->backward(right_tensor, &num_edmans, left_tensor);
    }
    double probability =
            backward_tensors.back()
                    .values[backward_tensors.back().strides[0] - 1];
    vector<Tensor>::iterator backward_tensor = backward_tensors.end();
    Tensor forward_tensor(tensor_shape.size(), &tensor_shape[0]);
    while (step != steps.end()) {
        backward_tensor--;
        (*step)->improve_fit(forward_tensor,
                             *backward_tensor,
                             *(backward_tensor - 1),
                             num_edmans,
                             probability,
                             fitter);
        (*step)->forward(forward_tensor, &num_edmans, &forward_tensor);
        step++;
    }
}

}  // namespace whatprot
