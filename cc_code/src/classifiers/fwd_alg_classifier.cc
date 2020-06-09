// Author: Matthew Beauregard Smith (UT Austin)
#include "fwd_alg_classifier.h"

#include <functional>

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

namespace fluoroseq {

FwdAlgClassifier::FwdAlgClassifier(
        int num_timesteps,
        int num_channels,
        const ErrorModel& error_model,
        const ApproximationModel& approximation_model,
        int num_dye_seqs,
        SourcedData<DyeSeq*, SourceWithCount<int>*>** dye_seqs)
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

}  // namespace fluoroseq
#include <iostream>
namespace fluoroseq {

ScoredClassification FwdAlgClassifier::classify(const Radiometry& radiometry) {
    Emission emission(radiometry, max_num_dyes, pdf, max_failed_edmans);
    int best_i = -1;
    double best_score = -1.0;
    double total_score = 0.0;
    // std::cout << "num_dye_seqs: " << num_dye_seqs << "\n";
    for (int i = 0; i < num_dye_seqs; i++) {
        // std::cout << "i: " << i << "\n";
        Initialization initialization;
        Summation summation(max_failed_edmans);
        double score = fwd_alg(tensors[i],
                               num_timesteps,
                               num_channels,
                               initialization,
                               emission,
                               *detach_transition,
                               *dud_transition,
                               *bleach_transition,
                               *edman_transitions[i],
                               summation);
        // std::cout << "score: " << score << "\n";
        total_score += score * dye_seqs[i]->source->count;
        if (score > best_score) {
            best_score = score;
            best_i = i;
        }
    }
    // std::cout << "best_i: " << best_i << "\n";
    auto a = dye_seqs[best_i];
    auto b = dye_seqs[best_i]->source;
    auto c = dye_seqs[best_i]->source->source;
    return ScoredClassification(dye_seqs[best_i]->source->source,
                                best_score,
                                total_score);
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
