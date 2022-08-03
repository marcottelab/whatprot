/******************************************************************************\
* Author: Matthew Beauregard Smith                                             *
* Affiliation: The University of Texas at Austin                               *
* Department: Oden Institute and Institute for Cellular and Molecular Biology  *
* PI: Edward Marcotte                                                          *
* Project: Protein Fluorosequencing                                            *
\******************************************************************************/

#ifndef WHATPROT_CLASSIFIERS_HMM_CLASSIFIER_H
#define WHATPROT_CLASSIFIERS_HMM_CLASSIFIER_H

// Standard C++ library headers:
#include <cmath>
#include <functional>
#include <vector>

// Local project headers:
#include "common/dye-seq.h"
#include "common/scored-classification.h"
#include "common/sourced-data.h"
#include "hmm/hmm/peptide-hmm.h"
#include "hmm/precomputations/dye-seq-precomputations.h"
#include "hmm/precomputations/radiometry-precomputations.h"
#include "hmm/precomputations/universal-precomputations.h"
#include "parameterization/model/sequencing-model.h"
#include "parameterization/settings/sequencing-settings.h"

namespace whatprot {

class HMMClassifier {
public:
    HMMClassifier(
            unsigned int num_timesteps,
            unsigned int num_channels,
            const SequencingModel& seq_model,
            const SequencingSettings& seq_settings,
            const std::vector<SourcedData<DyeSeq, SourceCount<int>>>& dye_seqs);
    ~HMMClassifier();
    ScoredClassification classify(const Radiometry& radiometry);
    ScoredClassification classify(const Radiometry& radiometry,
                                  const std::vector<int>& candidate_indices);
    std::vector<ScoredClassification> classify(
            const std::vector<Radiometry>& radiometries);
    std::vector<double> score(const Radiometry& radiometry);
    std::vector<std::vector<double>> score(
            const std::vector<Radiometry>& radiometries);

    template <typename I>
    ScoredClassification classify_helper(const Radiometry& radiometry,
                                         I indices) {
        RadiometryPrecomputations radiometry_precomputations(
                radiometry, seq_model, seq_settings, max_num_dyes);
        int best_i = -1;
        double best_score = -1.0;
        double total_score = 0.0;
        for (int i : indices) {
            PeptideHMM hmm(num_timesteps,
                           num_channels,
                           *dye_seq_precomputations_vec[i],
                           radiometry_precomputations,
                           universal_precomputations);
            double score = hmm.probability();
            total_score += score * dye_seqs[i].source.count;
            if (score > best_score) {
                best_score = score;
                best_i = i;
            }
        }
        ScoredClassification result(
                dye_seqs[best_i].source.source, best_score, total_score);
        // This next thing is a bit of a hack. Sometimes the candidates have a
        // total score of 0.0, which causes the adjusted score to be nan. This
        // can mess things up for us later. The best way to deal with it is to
        // just set the score to 0.0 when this happens. It might be better
        // though to find a way to avoid this situation.
        if (std::isnan(result.adjusted_score())) {
            result.score = 0.0;
            result.total = 1.0;
        }
        return result;
    }

    const SequencingModel& seq_model;
    const SequencingSettings& seq_settings;
    UniversalPrecomputations universal_precomputations;
    std::vector<DyeSeqPrecomputations*> dye_seq_precomputations_vec;
    const std::vector<SourcedData<DyeSeq, SourceCount<int>>>& dye_seqs;
    unsigned int num_timesteps;
    unsigned int num_channels;
    int max_num_dyes;
};

}  // namespace whatprot

#endif  // WHATPROT_CLASSIFIERS_HMM_CLASSIFIER_H
