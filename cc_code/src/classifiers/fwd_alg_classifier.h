// Author: Matthew Beauregard Smith (UT Austin)
#ifndef FLUOROSEQ_CLASSIFIERS_FWD_ALG_CLASSIFIER_H
#define FLUOROSEQ_CLASSIFIERS_FWD_ALG_CLASSIFIER_H

#include <functional>
#include <vector>

#include "common/approximation_model.h"
#include "common/dye_seq.h"
#include "common/error_model.h"
#include "common/scored_classification.h"
#include "common/sourced_data.h"
#include "fwd_alg/binomial_transition.h"
#include "fwd_alg/detach_transition.h"
#include "fwd_alg/edman_transition.h"
#include "fwd_alg/emission.h"
#include "fwd_alg/fwd_alg.h"
#include "fwd_alg/initialization.h"
#include "fwd_alg/summation.h"
#include "tensor/tensor.h"

namespace fluoroseq {

class FwdAlgClassifier {
public:
    FwdAlgClassifier(
            int num_timesteps,
            int num_channels,
            const ErrorModel& error_model,
            const ApproximationModel& approximation_model,
            const std::vector<SourcedData<DyeSeq, SourceCount<int>>>& dye_seqs);
    ~FwdAlgClassifier();
    ScoredClassification classify(const Radiometry& radiometry);
    ScoredClassification classify(const Radiometry& radiometry,
                                  const std::vector<int>& candidate_indices);
    std::vector<ScoredClassification> classify(
            const std::vector<Radiometry>& radiometries);

    template<typename I>
    ScoredClassification classify_helper(const Radiometry& radiometry,
                                         I indices) {
        Emission emission(radiometry, max_num_dyes, pdf, max_failed_edmans);
        int best_i = -1;
        double best_score = -1.0;
        double total_score = 0.0;
        int i = 0;
        for (int i : indices) {
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
            total_score += score * dye_seqs[i].source.count;
            if (score > best_score) {
                best_score = score;
                best_i = i;
            }
        }
        return ScoredClassification(dye_seqs[best_i].source.source,
                                    best_score,
                                    total_score);

    }

    DetachTransition* detach_transition;
    BinomialTransition* dud_transition;
    BinomialTransition* bleach_transition;
    std::function<double (double, int)> pdf;
    const std::vector<SourcedData<DyeSeq, SourceCount<int>>>& dye_seqs;
    EdmanTransition** edman_transitions;
    Tensor** tensors;
    int num_dye_seqs;
    int num_timesteps;
    int num_channels;
    int max_num_dyes;
    int max_failed_edmans;
};

}  // namespace fluoroseq

#endif  // FLUOROSEQ_CLASSIFIERS_FWD_ALG_CLASSIFIER_H
