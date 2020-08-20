// Author: Matthew Beauregard Smith (UT Austin)
#ifndef FLUOROSEQ_CLASSIFIERS_HYBRID_CLASSIFIER_H
#define FLUOROSEQ_CLASSIFIERS_HYBRID_CLASSIFIER_H

#include <unordered_map>
#include <vector>

#include "classifiers/fwd_alg_classifier.h"
#include "classifiers/kwann_classifier.h"
#include "common/approximation_model.h"
#include "common/dye_seq.h"
#include "common/dye_track.h"
#include "common/error_model.h"
#include "common/radiometry.h"
#include "common/scored_classification.h"
#include "common/sourced_data.h"

namespace fluoroseq {

class HybridClassifier {
public:
    HybridClassifier(
            int num_timesteps,
            int num_channels,
            const ErrorModel& error_model,
            const ApproximationModel& approximation_model,
            int k,
            const std::vector<
                    SourcedData<DyeTrack,
                                SourceCountHitsList<int>>>& dye_tracks,
            int h,
            const std::vector<SourcedData<DyeSeq, SourceCount<int>>>& dye_seqs);
    ~HybridClassifier();
    ScoredClassification classify(const Radiometry& radiometry);
    std::vector<ScoredClassification> classify(
            const std::vector<Radiometry>& radiometries);

    FwdAlgClassifier* fwd_alg_classifier;
    KWANNClassifier* kwann_classifier;
    std::unordered_map<int, int> id_index_map;
    std::unordered_map<int, int> id_count_map;
    int h;
};

}  // namespace fluoroseq

#endif  // FLUOROSEQ_CLASSIFIERS_HYBRID_CLASSIFIER_H