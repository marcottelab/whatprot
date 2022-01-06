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

namespace {
using std::isnan;
using std::vector;
}  // namespace

HybridClassifier::HybridClassifier(
        unsigned int num_timesteps,
        unsigned int num_channels,
        const SequencingModel& seq_model,
        const SequencingSettings& seq_settings,
        int k,
        double sig,
        vector<SourcedData<DyeTrack, SourceCountHitsList<int>>>* dye_tracks,
        int h,
        const vector<SourcedData<DyeSeq, SourceCount<int>>>& dye_seqs)
        : hmm_classifier(
                num_timesteps, num_channels, seq_model, seq_settings, dye_seqs),
          nn_classifier(
                  num_timesteps, num_channels, seq_model, k, sig, dye_tracks),
          h(h) {
    for (unsigned int i = 0; i < dye_seqs.size(); i++) {
        id_index_map[dye_seqs[i].source.source] = i;
        id_count_map[dye_seqs[i].source.source] = dye_seqs[i].source.count;
    }
}

ScoredClassification HybridClassifier::classify(const Radiometry& radiometry) {
    vector<ScoredClassification> candidates;
    candidates = nn_classifier.classify(radiometry, h);
    double subfraction = 0.0;
    vector<int> candidate_indices;
    candidate_indices.reserve(candidates.size());
    for (ScoredClassification& candidate : candidates) {
        subfraction +=
                candidate.adjusted_score() * (double)id_count_map[candidate.id];
        candidate_indices.push_back(id_index_map[candidate.id]);
    }
    ScoredClassification result;
    result = hmm_classifier.classify(radiometry, candidate_indices);
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
    results.resize(radiometries.size());
#pragma omp parallel for
    for (unsigned int i = 0; i < radiometries.size(); i++) {
        results[i] = classify(radiometries[i]);
    }
    return results;
}

}  // namespace whatprot
