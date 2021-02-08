/******************************************************************************\
* Author: Matthew Beauregard Smith                                             *
* Affiliation: The University of Texas at Austin                               *
* Department: Oden Institute and Institute for Cellular and Molecular Biology  *
* PI: Edward Marcotte                                                          *
* Project: Protein Fluorosequencing                                            *
\******************************************************************************/

// Defining symbols from header:
#include "fwd-alg-classifier.h"

// Standard C++ library headers:
#include <functional>
#include <utility>
#include <vector>

// Local project headers:
#include "common/dye-seq.h"
#include "common/dye-track.h"
#include "common/error-model.h"
#include "common/scored-classification.h"
#include "fwd-alg/binomial-transition.h"
#include "fwd-alg/detach-transition.h"
#include "fwd-alg/edman-transition.h"
#include "fwd-alg/emission.h"
#include "fwd-alg/fwd-alg.h"
#include "fwd-alg/initialization.h"
#include "fwd-alg/summation.h"
#include "util/range.h"

namespace whatprot {

namespace {
using std::function;
using std::move;
using std::vector;
}  // namespace

FwdAlgClassifier::FwdAlgClassifier(
        int num_timesteps,
        int num_channels,
        const ErrorModel& error_model,
        const vector<SourcedData<DyeSeq, SourceCount<int>>>& dye_seqs)
        : num_timesteps(num_timesteps),
          num_channels(num_channels),
          num_dye_seqs(dye_seqs.size()),
          dye_seqs(dye_seqs),
          detach_transition(error_model.p_detach),
          dud_transition(error_model.p_bleach),
          bleach_transition(error_model.p_bleach) {
    edman_transitions.reserve(num_dye_seqs);
    tensor_shapes.resize(num_dye_seqs);
    max_num_dyes = 0;
    for (int i = 0; i < num_dye_seqs; i++) {
        DyeTrack dye_track =
                DyeTrack(num_timesteps, num_channels, dye_seqs[i].value);
        for (int c = 0; c < num_channels; c++) {
            if (dye_track(0, c) > max_num_dyes) {
                max_num_dyes = dye_track(0, c);
            }
        }
        edman_transitions.push_back(move(EdmanTransition(
                error_model.p_edman_failure, dye_seqs[i].value, dye_track)));
        tensor_shapes[i].resize(1 + num_channels);
        tensor_shapes[i][0] = num_timesteps + 1;
        for (int c = 0; c < num_channels; c++) {
            // This next line of code is a little confusing.
            //   * The first dimension of the tensor shape is always the
            //     timestep, so we need to add one to the channel to index to
            //     the correct dimension.
            //   * The zeroth step for the dye track gives us the maximum number
            //     of dyes possible, but the tensor shape for that channel needs
            //     to be one bigger than that to handle all values inclusively
            //     from 0 to the number of dyes.
            tensor_shapes[i][1 + c] = 1 + dye_track(0, c);
        }
    }
    dud_transition.reserve(max_num_dyes);
    bleach_transition.reserve(max_num_dyes);
    pdf = error_model.pdf();
}

ScoredClassification FwdAlgClassifier::classify(const Radiometry& radiometry) {
    return classify_helper<Range>(radiometry, Range(num_dye_seqs));
}

ScoredClassification FwdAlgClassifier::classify(
        const Radiometry& radiometry, const vector<int>& candidate_indices) {
    return classify_helper<const vector<int>&>(radiometry, candidate_indices);
}

vector<ScoredClassification> FwdAlgClassifier::classify(
        const vector<Radiometry>& radiometries) {
    vector<ScoredClassification> results;
    results.resize(radiometries.size());
    #pragma omp parallel for
    for (int i = 0; i < radiometries.size(); i++) {
        results[i] = classify(radiometries[i]);
    }
    return results;
}

}  // namespace whatprot
