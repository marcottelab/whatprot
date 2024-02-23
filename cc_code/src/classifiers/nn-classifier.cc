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
#include "common/radiometry.h"
#include "common/scored-classification.h"
#include "common/sourced-data.h"
#include "parameterization/model/sequencing-model.h"

namespace {
using std::exp;
using std::function;
using std::greater;  // defined in <functional>
using std::isnan;
using std::move;
using std::priority_queue;
using std::sqrt;
using std::unordered_map;
using std::vector;
using whatprot::KDTEntry;  // in namespace std for swap
}  // namespace

namespace whatprot {

KDTEntry::KDTEntry(const SequencingModel& seq_model,
                   SourcedData<DyeTrack, SourceCountHitsList<int>>&& dye_track)
        : seq_model(&seq_model), dye_track(move(dye_track)) {
    hits = this->dye_track.source.total_hits();
}

KDTEntry::KDTEntry(KDTEntry&& other)
        : seq_model(other.seq_model),
          dye_track(move(other.dye_track)),
          hits(other.hits) {}

KDTEntry& KDTEntry::operator=(KDTEntry&& other) {
    seq_model = other.seq_model;
    dye_track = move(other.dye_track);
    hits = other.hits;
    return *this;
}

double KDTEntry::operator[](int i) const {
    unsigned int num_channels = dye_track.value.num_channels;
    // The dye_track counts are served by a single index which flattens the
    // indexing scheme. We extract c (the channel) and t (the timestep).
    unsigned int c = i % num_channels;
    unsigned int t = i / num_channels;
    // To compute the mu value accounting for cross-dye interactions, we need
    // counts for all channels at the timestep. This is a bit ugly but we can
    // get a pointer to the value for channel 0 at the desired timestep t (by
    // indexing at t * num_channels), and adjusted_mu() will index the various
    // channels as needed.
    const short* t_counts = &dye_track.value.counts[t * num_channels];
    // Now we can use the channel_model for our channel of interest (c) and
    // retrieve the adjusted_mu value based on t_counts.
    return seq_model->channel_models[c]->adjusted_mu(t_counts);
}

}  // namespace whatprot

namespace std {

void swap(KDTEntry& e1, KDTEntry& e2) {
    KDTEntry temp(move(e1));
    e1 = move(e2);
    e2 = move(temp);
}

}  // namespace std

namespace whatprot {

KDTQuery::KDTQuery(const Radiometry& rad) : rad(rad) {}

double KDTQuery::operator[](int i) const {
    return rad.intensities[i];
}

NNClassifier::NNClassifier(
        unsigned int num_timesteps,
        unsigned int num_channels,
        const SequencingModel& seq_model,
        int k,
        double sig,
        vector<SourcedData<DyeTrack, SourceCountHitsList<int>>>* dye_tracks)
        : num_train(dye_tracks->size()),
          num_timesteps(num_timesteps),
          num_channels(num_channels),
          k(k),
          two_sig_sq(2.0 * sig * sig) {
    vector<KDTEntry> kdt_entries;
    kdt_entries.reserve(num_train);
    for (int i = 0; i < num_train; i++) {
        KDTEntry kdt_convert(seq_model, move((*dye_tracks)[i]));
        kdt_entries.push_back(move(kdt_convert));
    }
    kd_tree = new KDTree<KDTEntry, KDTQuery>(k,
                                             num_timesteps * num_channels,  // d
                                             move(kdt_entries));
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
    for (unsigned int i = 0; i < k_nearest.size(); i++) {
        const SourcedData<DyeTrack, SourceCountHitsList<int>>& dye_track =
                k_nearest[i]->dye_track;
        double dist_sq = dists_sq[i];
        // For computing a gaussian kernel.
        //   * The normalization factor, 1/(sig*2*PI), is ignored here,
        //     because it is a constant factor, so all weights should be
        //     affected equally.
        //   * We use the dist_sq from the KDTree. This works because a guassian
        //     kernel is radially symmetric. It is also far more efficient to
        //     compute it this way, which is why we do it.
        double weight = exp(-dist_sq / two_sig_sq);
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
        const Radiometry& radiometry, unsigned int h) {
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
    results.resize(radiometries.size());
#pragma omp parallel for schedule(dynamic, 1)
    for (unsigned int i = 0; i < radiometries.size(); i++) {
        results[i] = classify(radiometries[i]);
    }
    return results;
}

}  // namespace whatprot
