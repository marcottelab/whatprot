/******************************************************************************\
* Author: Matthew Beauregard Smith                                             *
* Affiliation: The University of Texas at Austin                               *
* Department: Oden Institute and Institute for Cellular and Molecular Biology  *
* PI: Edward Marcotte                                                          *
* Project: Protein Fluorosequencing                                            *
\******************************************************************************/

#ifndef FLUOROSEQ_CLASSIFIERS_HMM_CLASSIFIER_H
#define FLUOROSEQ_CLASSIFIERS_HMM_CLASSIFIER_H

// Standard C++ library headers:
#include <functional>
#include <vector>

// Local project headers:
#include "common/dye-seq.h"
#include "common/error-model.h"
#include "common/scored-classification.h"
#include "common/sourced-data.h"
#include "hmm/hmm.h"
#include "hmm/precomputations/dye-seq-precomputations.h"
#include "hmm/precomputations/radiometry-precomputations.h"
#include "hmm/precomputations/universal-precomputations.h"

namespace fluoroseq {

class HMMClassifier {
public:
    HMMClassifier(
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
        RadiometryPrecomputations radiometry_precomputations(
                radiometry, error_model, max_num_dyes);
        int best_i = -1;
        double best_score = -1.0;
        double total_score = 0.0;
        int i = 0;
        for (int i : indices) {
            HMM hmm(num_timesteps,
                    num_channels,
                    dye_seq_precomputations_vec[i],
                    radiometry_precomputations,
                    universal_precomputations);
            double score = hmm.probability();
            total_score += score * dye_seqs[i].source.count;
            if (score > best_score) {
                best_score = score;
                best_i = i;
            }
        }
        return ScoredClassification(
                dye_seqs[best_i].source.source, best_score, total_score);
    }

    const ErrorModel& error_model;
    UniversalPrecomputations universal_precomputations;
    std::vector<DyeSeqPrecomputations> dye_seq_precomputations_vec;
    const std::vector<SourcedData<DyeSeq, SourceCount<int>>>& dye_seqs;
    int num_timesteps;
    int num_channels;
    int max_num_dyes;
};

}  // namespace fluoroseq

#endif  // FLUOROSEQ_CLASSIFIERS_HMM_CLASSIFIER_H
