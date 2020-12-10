/******************************************************************************\
* Author: Matthew Beauregard Smith                                             *
* Affiliation: The University of Texas at Austin                               *
* Department: Oden Institute and Institute for Cellular and Molecular Biology  *
* PI: Edward Marcotte                                                          *
* Project: Protein Fluorosequencing                                            *
\******************************************************************************/

// Defining symbols from header:
#include "nn-classifier.h"

// Standard C++ library headers:
#include <cmath>
#include <functional>
#include <queue>
#include <unordered_map>
#include <utility>
#include <vector>

// Local project headers:
#include "common/error-model.h"
#include "common/radiometry.h"
#include "common/scored-classification.h"
#include "common/sourced-data.h"

namespace {
using fluoroseq::KDTEntry;  // in namespace std for swap
using std::exp;
using std::function;
using std::greater;  // defined in <functional>
using std::isnan;
using std::move;
using std::priority_queue;
using std::sqrt;
using std::unordered_map;
using std::vector;
double PI = 3.141592653589793238;
}  // namespace

namespace fluoroseq {

KDTEntry::KDTEntry(SourcedData<DyeTrack, SourceCountHitsList<int>>&& dye_track)
        : dye_track(move(dye_track)) {}

KDTEntry::KDTEntry(KDTEntry&& other) : dye_track(move(other.dye_track)) {}

KDTEntry& KDTEntry::operator=(KDTEntry&& other) {
    dye_track = move(other.dye_track);
    return *this;
}

double KDTEntry::operator[](int i) const {
    return (double)dye_track.value.counts[i];
}

}  // namespace fluoroseq

namespace std {

void swap(KDTEntry& e1, KDTEntry& e2) {
    KDTEntry temp(move(e1));
    e1 = move(e2);
    e2 = move(temp);
}

}  // namespace std

namespace fluoroseq {

KDTQuery::KDTQuery(const Radiometry& rad) : rad(rad) {}

double KDTQuery::operator[](int i) const {
    return rad.intensities[i];
}

NNClassifier::NNClassifier(
        int num_timesteps,
        int num_channels,
        const ErrorModel& error_model,
        int k,
        double sigma,
        vector<SourcedData<DyeTrack, SourceCountHitsList<int>>>* dye_tracks)
        : num_timesteps(num_timesteps),
          num_channels(num_channels),
          k(k),
          num_train(dye_tracks->size()) {
    int stride = num_timesteps * num_channels;
    vector<KDTEntry> kdt_entries;
    kdt_entries.reserve(num_train);
    for (int i = 0; i < num_train; i++) {
        KDTEntry kdt_convert(move((*dye_tracks)[i]));
        kdt_entries.push_back(move(kdt_convert));
    }
    kd_tree = new KDTree<KDTEntry, KDTQuery>(k,
                                             num_timesteps * num_channels,  // d
                                             move(kdt_entries));
    double scale = error_model.mu;
    double sig = sigma;
    double multiplier = 1.0 / (sigma * sqrt(2.0 * PI));
    kernel = [scale, sig, multiplier](double observed, int state) -> double {
        double unit_obs = observed / scale;
        double offset = unit_obs - (double)state;
        return multiplier * exp(-(offset * offset) / (2.0 * sig * sig));
    };
}

NNClassifier::~NNClassifier() {
    delete kd_tree;
}

double NNClassifier::classify_helper(const Radiometry& radiometry,
                                     unordered_map<int, double>* id_score_map) {
    KDTQuery query(radiometry);
    vector<KDTEntry*> k_nearest;
    vector<double> dists_sq;
    kd_tree->search(query, &k_nearest, &dists_sq);
    double total_score = 0.0;
    for (int i = 0; i < k; i++) {
        const SourcedData<DyeTrack, SourceCountHitsList<int>>& dye_track =
                k_nearest[i]->dye_track;
        double weight = 1.0;
        for (int j = 0; j < num_timesteps * num_channels; j++) {
            double offset = radiometry.intensities[j]
                            - (double)dye_track.value.counts[j];
            weight *= kernel(radiometry.intensities[j],
                             dye_track.value.counts[j]);
        }
        for (int j = 0; j < dye_track.source.num_sources; j++) {
            int id = dye_track.source.sources[j]->source;
            double count = (double)dye_track.source.sources[j]->count;
            double hits = (double)dye_track.source.sources[j]->hits;
            total_score += weight * hits;
            (*id_score_map)[id] += weight * hits / count;
        }
    }
    return total_score;
}

ScoredClassification NNClassifier::classify(const Radiometry& radiometry) {
    unordered_map<int, double> id_score_map;
    double total_score = classify_helper(radiometry, &id_score_map);
    int best_id = -1;
    double best_score = -1.0;
    for (const auto& id_and_score : id_score_map) {
        double id = id_and_score.first;
        double score = id_and_score.second;
        if (score > best_score) {
            best_id = id;
            best_score = score;
        }
    }
    ScoredClassification result(best_id, best_score, total_score);
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

vector<ScoredClassification> NNClassifier::classify(
        const Radiometry& radiometry, int h) {
    unordered_map<int, double> id_score_map;
    double total_score = classify_helper(radiometry, &id_score_map);
    priority_queue<ScoredClassification,
                   vector<ScoredClassification>,
                   greater<ScoredClassification>>
            pq;
    for (const auto& id_and_score : id_score_map) {
        double id = id_and_score.first;
        double score = id_and_score.second;
        if (pq.size() < h) {
            pq.push(ScoredClassification(id, score, total_score));
        } else if (score > pq.top().score) {
            pq.push(ScoredClassification(id, score, total_score));
            pq.pop();
        }
    }
    vector<ScoredClassification> results;
    results.reserve(pq.size());
    while (!pq.empty()) {
        results.push_back(pq.top());
        pq.pop();
    }
    return results;
}

vector<ScoredClassification> NNClassifier::classify(
        const vector<Radiometry>& radiometries) {
    vector<ScoredClassification> results;
    results.reserve(radiometries.size());
    for (int i = 0; i < radiometries.size(); i++) {
        results.push_back(classify(radiometries[i]));
    }
    return results;
}

}  // namespace fluoroseq
