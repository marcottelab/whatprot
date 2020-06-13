// Author: Matthew Beauregard Smith (UT Austin)
#include "hybrid_classifier.h"

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

namespace {
using std::vector;
}

HybridClassifier::HybridClassifier(
        int num_timesteps,
        int num_channels,
        const ErrorModel& error_model,
        const ApproximationModel& approximation_model,
        int k,
        int num_train,
        SourcedData<DyeTrack*, SourceCountMap<int>*>** dye_tracks,
        int h,
        int num_dye_seqs,
        SourcedData<DyeSeq*, SourceWithCount<int>*>** dye_seqs) : h(h) {
    kwann_classifier = new KWANNClassifier(num_timesteps,
                                           num_channels,
                                           error_model.pdf(),
                                           k,
                                           num_train,
                                           dye_tracks);
    fwd_alg_classifier = new FwdAlgClassifier(num_timesteps,
                                              num_channels,
                                              error_model,
                                              approximation_model,
                                              num_dye_seqs,
                                              dye_seqs);
    for (int i = 0; i < num_dye_seqs; i++) {
        id_index_map[dye_seqs[i]->source->source] = i;
        id_count_map[dye_seqs[i]->source->source] = dye_seqs[i]->source->count;
    }
}

HybridClassifier::~HybridClassifier() {
    delete kwann_classifier;
    delete fwd_alg_classifier;
}

}  // namespace_fluoroseq
#include <iostream>
namespace fluoroseq {

ScoredClassification HybridClassifier::classify(const Radiometry& radiometry) {
    vector<ScoredClassification> candidates;
    candidates = kwann_classifier->classify(radiometry, h);
    std::cout << "candidates.size(): " << candidates.size() << "\n";
    double total = candidates[0].total;
    double subtotal = 0.0;
    vector<int> candidate_indices;
    candidate_indices.reserve(candidates.size());
    for (ScoredClassification& candidate : candidates) {
        subtotal += candidate.score * id_count_map[candidate.id];
        candidate_indices.push_back(id_index_map[candidate.id]);
    }
    double total_correction_ratio = total / subtotal;
    ScoredClassification result;
    result = fwd_alg_classifier->classify(radiometry, candidate_indices);
    if (result.id == -1) {
        result = candidates.back();
    } else {
        result.total *= total_correction_ratio;
    }
    return result;
}

ScoredClassification* HybridClassifier::classify(int num_radiometries, 
                                                 Radiometry** radiometries) {
    ScoredClassification* results = new ScoredClassification[num_radiometries];
    for (int i = 0; i < num_radiometries; i++) {
        results[i] = classify(*radiometries[i]);
    }
    return results;
}

}  // namespace fluoroseq
