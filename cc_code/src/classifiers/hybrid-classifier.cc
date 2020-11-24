/******************************************************************************\
* Author: Matthew Beauregard Smith                                             *
* Affiliation: The University of Texas at Austin                               *
* Department: Oden Institute and Institute for Cellular and Molecular Biology  *
* PI: Edward Marcotte                                                          *
* Project: Protein Fluorosequencing                                            *
\******************************************************************************/

// Defining symbols from header:
#include "hybrid-classifier.h"

// Standard C++ library headers:
#include <cmath>
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

namespace {
using std::isnan;
using std::vector;
}  // namespace

HybridClassifier::HybridClassifier(
        int num_timesteps,
        int num_channels,
        const ErrorModel& error_model,
        int k,
        double sigma,
        const vector<SourcedData<DyeTrack, SourceCountHitsList<int>>>&
                dye_tracks,
        int h,
        const vector<SourcedData<DyeSeq, SourceCount<int>>>& dye_seqs)
        : h(h),
          kwann_classifier(num_timesteps,
                           num_channels,
                           error_model,
                           k,
                           sigma,
                           dye_tracks),
          fwd_alg_classifier(
                  num_timesteps, num_channels, error_model, dye_seqs) {
    for (int i = 0; i < dye_seqs.size(); i++) {
        id_index_map[dye_seqs[i].source.source] = i;
        id_count_map[dye_seqs[i].source.source] = dye_seqs[i].source.count;
    }
}

ScoredClassification HybridClassifier::classify(const Radiometry& radiometry) {
    vector<ScoredClassification> candidates;
    candidates = kwann_classifier.classify(radiometry, h);
    double total = candidates[0].total;
    double subfraction = 0.0;
    vector<int> candidate_indices;
    candidate_indices.reserve(candidates.size());
    for (ScoredClassification& candidate : candidates) {
        subfraction +=
                candidate.adjusted_score() * (double)id_count_map[candidate.id];
        candidate_indices.push_back(id_index_map[candidate.id]);
    }
    ScoredClassification result;
    result = fwd_alg_classifier.classify(radiometry, candidate_indices);
    if (result.id == -1) {
        result = candidates.back();
    } else {
        result.score *= subfraction;
    }
    // This next thing is a bit of a hack. Sometimes the candidates have a total
    // score of 0.0, which causes the adjusted score to be nan. This can mess
    // things up for us later. The best way to deal with it is to just set the
    // score to 0.0 when this happens. It might be better though to find a way
    // to avoid this situation.
    if (isnan(result.adjusted_score())) {
        result.score = 0.0;
        result.total = 1.0;
    }
    return result;
}

vector<ScoredClassification> HybridClassifier::classify(
        const vector<Radiometry>& radiometries) {
    vector<ScoredClassification> results;
    results.reserve(radiometries.size());
    for (int i = 0; i < radiometries.size(); i++) {
        results.push_back(classify(radiometries[i]));
    }
    return results;
}

}  // namespace fluoroseq