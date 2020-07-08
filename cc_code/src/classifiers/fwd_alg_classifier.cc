// Author: Matthew Beauregard Smith (UT Austin)
#include "fwd_alg_classifier.h"

#include <functional>
#include <vector>

#include "common/dye_seq.h"
#include "common/dye_track.h"
#include "common/error_model.h"
#include "common/scored_classification.h"
#include "fwd_alg/binomial_transition.h"
#include "fwd_alg/detach_transition.h"
#include "fwd_alg/edman_transition.h"
#include "fwd_alg/emission.h"
#include "fwd_alg/fwd_alg.h"
#include "fwd_alg/initialization.h"
#include "fwd_alg/summation.h"
#include "util/range.h"

namespace fluoroseq {

namespace {
using std::function;
using std::vector;
}

FwdAlgClassifier::FwdAlgClassifier(
        int num_timesteps,
        int num_channels,
        const ErrorModel& error_model,
        const ApproximationModel& approximation_model,
        int num_dye_seqs,
        SourcedData<DyeSeq*, SourceCount<int>*>** dye_seqs)
        : num_timesteps(num_timesteps),
          num_channels(num_channels),
          num_dye_seqs(num_dye_seqs),
          dye_seqs(dye_seqs),
          max_failed_edmans(approximation_model.max_failed_edmans) {
    detach_transition = new DetachTransition(error_model.p_detach,
                                             max_failed_edmans);
    edman_transitions = new EdmanTransition*[num_dye_seqs]();
    tensors = new Tensor*[num_dye_seqs]();
    max_num_dyes = 0;
    for (int i = 0; i < num_dye_seqs; i++) {
        DyeTrack dye_track = DyeTrack(num_timesteps,
                                      num_channels,
                                      *dye_seqs[i]->value);
        for (int c = 0; c < num_channels; c++) {
            if (dye_track(0, c) > max_num_dyes) {
                max_num_dyes = dye_track(0, c);
            }
        }
        edman_transitions[i] = new EdmanTransition(error_model.p_edman_failure,
                                                   *dye_seqs[i]->value,
                                                   dye_track,
                                                   max_failed_edmans);
        int* tensor_shape = new int[1 + num_channels];
        tensor_shape[0] = num_timesteps + 1;
        for (int c = 0; c < num_channels; c++) {
            int num_dyes = dye_track(0, c);
            tensor_shape[1 + c] = num_dyes + 1;
        }
        tensors[i] = new Tensor(1 + num_channels, tensor_shape);
        delete[] tensor_shape;
    }
    dud_transition = new BinomialTransition(max_num_dyes,
                                            error_model.p_dud,
                                            max_failed_edmans);
    bleach_transition = new BinomialTransition(max_num_dyes,
                                               error_model.p_bleach,
                                               max_failed_edmans);
    pdf = error_model.pdf();
}

FwdAlgClassifier::~FwdAlgClassifier() {
    delete detach_transition;
    delete dud_transition;
    delete bleach_transition;
    for (int i = 0; i < num_dye_seqs; i++) {
        delete edman_transitions[i];
    }
    delete[] edman_transitions;
    for (int i = 0; i < num_dye_seqs; i++) {
        delete tensors[i];
    }
    delete[] tensors;
}

ScoredClassification FwdAlgClassifier::classify(const Radiometry& radiometry) {
    return classify_helper<Range>(radiometry, Range(num_dye_seqs));
}

ScoredClassification FwdAlgClassifier::classify(
        const Radiometry& radiometry,
        const vector<int>& candidate_indices) {
    return classify_helper<const vector<int>&>(radiometry, candidate_indices);
}

ScoredClassification* FwdAlgClassifier::classify(int num_radiometries,
                                                 Radiometry** radiometries) {
    ScoredClassification* result = new ScoredClassification[num_radiometries];
    for (int i = 0; i < num_radiometries; i++) {
        result[i] = classify(*radiometries[i]);
    }
    return result;
}

}  // namespace fluoroseq
