/******************************************************************************\
* Author: Matthew Beauregard Smith                                             *
* Affiliation: The University of Texas at Austin                               *
* Department: Oden Institute and Institute for Cellular and Molecular Biology  *
* PI: Edward Marcotte                                                          *
* Project: Protein Fluorosequencing                                            *
\******************************************************************************/

#ifndef FLUOROSEQ_CLASSIFIERS_HYBRID_CLASSIFIER_H
#define FLUOROSEQ_CLASSIFIERS_HYBRID_CLASSIFIER_H

// Standard C++ library headers:
#include <unordered_map>
#include <vector>

// Local project headers:
#include "classifiers/fwd-alg-classifier.h"
#include "classifiers/kwann-classifier.h"
#include "common/dye-seq.h"
#include "common/dye-track.h"
#include "common/error-model.h"
#include "common/radiometry.h"
#include "common/scored-classification.h"
#include "common/sourced-data.h"

namespace fluoroseq {

class HybridClassifier {
public:
    HybridClassifier(
            int num_timesteps,
            int num_channels,
            const ErrorModel& error_model,
            int k,
            double sigma,
            const std::vector<SourcedData<DyeTrack, SourceCountHitsList<int>>>&
                    dye_tracks,
            int h,
            const std::vector<SourcedData<DyeSeq, SourceCount<int>>>& dye_seqs);
    ScoredClassification classify(const Radiometry& radiometry);
    std::vector<ScoredClassification> classify(
            const std::vector<Radiometry>& radiometries);

    FwdAlgClassifier fwd_alg_classifier;
    KWANNClassifier kwann_classifier;
    std::unordered_map<int, int> id_index_map;
    std::unordered_map<int, int> id_count_map;
    int h;
};

}  // namespace fluoroseq

#endif  // FLUOROSEQ_CLASSIFIERS_HYBRID_CLASSIFIER_H