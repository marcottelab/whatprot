/******************************************************************************\
* Author: Matthew Beauregard Smith                                             *
* Affiliation: The University of Texas at Austin                               *
* Department: Oden Institute and Institute for Cellular and Molecular Biology  *
* PI: Edward Marcotte                                                          *
* Project: Protein Fluorosequencing                                            *
\******************************************************************************/

#ifndef WHATPROT_CLASSIFIERS_HYBRID_CLASSIFIER_H
#define WHATPROT_CLASSIFIERS_HYBRID_CLASSIFIER_H

// Standard C++ library headers:
#include <unordered_map>
#include <vector>

// Local project headers:
#include "classifiers/hmm-classifier.h"
#include "classifiers/nn-classifier.h"
#include "common/dye-seq.h"
#include "common/dye-track.h"
#include "common/radiometry.h"
#include "common/scored-classification.h"
#include "common/sourced-data.h"
#include "parameterization/model/sequencing-model.h"
#include "parameterization/settings/sequencing-settings.h"

namespace whatprot {

class HybridClassifier {
public:
    HybridClassifier(
            unsigned int num_timesteps,
            unsigned int num_channels,
            const SequencingModel& seq_model,
            const SequencingSettings& seq_settings,
            int k,
            double sigma,
            std::vector<SourcedData<DyeTrack, SourceCountHitsList<int>>>*
                    dye_tracks,
            int h,
            const std::vector<SourcedData<DyeSeq, SourceCount<int>>>& dye_seqs);
    ScoredClassification classify(const Radiometry& radiometry);
    std::vector<ScoredClassification> classify(
            const std::vector<Radiometry>& radiometries);

    HMMClassifier hmm_classifier;
    NNClassifier nn_classifier;
    std::unordered_map<int, int> id_index_map;
    std::unordered_map<int, int> id_count_map;
    int h;
};

}  // namespace whatprot

#endif  // WHATPROT_CLASSIFIERS_HYBRID_CLASSIFIER_H