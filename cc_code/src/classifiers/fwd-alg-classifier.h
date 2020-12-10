/******************************************************************************\
* Author: Matthew Beauregard Smith                                             *
* Affiliation: The University of Texas at Austin                               *
* Department: Oden Institute and Institute for Cellular and Molecular Biology  *
* PI: Edward Marcotte                                                          *
* Project: Protein Fluorosequencing                                            *
\******************************************************************************/

#ifndef FLUOROSEQ_CLASSIFIERS_FWD_ALG_CLASSIFIER_H
#define FLUOROSEQ_CLASSIFIERS_FWD_ALG_CLASSIFIER_H

// Standard C++ library headers:
#include <functional>
#include <vector>

// Local project headers:
#include "common/dye-seq.h"
#include "common/error-model.h"
#include "common/scored-classification.h"
#include "common/sourced-data.h"
#include "fwd-alg/binomial-transition.h"
#include "fwd-alg/detach-transition.h"
#include "fwd-alg/edman-transition.h"
#include "fwd-alg/emission.h"
#include "fwd-alg/fwd-alg.h"
#include "fwd-alg/initialization.h"
#include "fwd-alg/summation.h"
#include "tensor/tensor.h"

namespace fluoroseq {

class FwdAlgClassifier {
public:
    FwdAlgClassifier(
            int num_timesteps,
            int num_channels,
            const ErrorModel& error_model,
            const std::vector<SourcedData<DyeSeq, SourceCount<int>>>& dye_seqs);
    ScoredClassification classify(const Radiometry& radiometry);
    ScoredClassification classify(const Radiometry& radiometry,
                                  const std::vector<int>& candidate_indices);
    std::vector<ScoredClassification> classify(
            const std::vector<Radiometry>& radiometries);

    template <typename I>
    ScoredClassification classify_helper(const Radiometry& radiometry,
                                         I indices) {
        Emission emission(radiometry, max_num_dyes, pdf);
        int best_i = -1;
        double best_score = -1.0;
        double total_score = 0.0;
        int i = 0;
        for (int i : indices) {
            Initialization initialization;
            Summation summation;
            double score = fwd_alg(&tensors[i],
                                   num_timesteps,
                                   num_channels,
                                   initialization,
                                   emission,
                                   detach_transition,
                                   dud_transition,
                                   bleach_transition,
                                   edman_transitions[i],
                                   summation);
            total_score += score * dye_seqs[i].source.count;
            if (score > best_score) {
                best_score = score;
                best_i = i;
            }
        }
        return ScoredClassification(
                dye_seqs[best_i].source.source, best_score, total_score);
    }

    DetachTransition detach_transition;
    BinomialTransition dud_transition;
    BinomialTransition bleach_transition;
    std::function<double(double, int)> pdf;
    const std::vector<SourcedData<DyeSeq, SourceCount<int>>>& dye_seqs;
    std::vector<EdmanTransition> edman_transitions;
    std::vector<Tensor> tensors;
    int num_dye_seqs;
    int num_timesteps;
    int num_channels;
    int max_num_dyes;
};

}  // namespace fluoroseq

#endif  // FLUOROSEQ_CLASSIFIERS_FWD_ALG_CLASSIFIER_H
