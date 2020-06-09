// Author: Matthew Beauregard Smith (UT Austin)
#include "kwann_classifier.h"

#include <functional>
#include <unordered_map>

#include "common/radiometry.h"
#include "common/scored_classification.h"
#include "common/sourced_data.h"
#include "flann/flann.hpp"

namespace {
using flann::AutotunedIndexParams;
using flann::FLANN_CHECKS_AUTOTUNED;
using flann::Index;
using flann::L2;
using flann::Matrix;
using flann::SearchParams;
using std::function;
using std::unordered_map;
}  // namespace

namespace fluoroseq {

KWANNClassifier::KWANNClassifier(int num_timesteps,
                                 int num_channels,
                                 function<double (double, int)> pdf,
                                 int k,
                                 int num_train,
                                 SourcedData<DyeTrack*,
                                             SourceCountMap<int>*>** dye_tracks)
        : num_timesteps(num_timesteps),
          num_channels(num_channels),
          pdf(pdf),
          k(k),
          num_train(num_train),
          dye_tracks(dye_tracks) {
    int stride = num_timesteps * num_channels;
    double* raw_dataset = new double[num_train * stride];
    for (int i = 0; i < num_train; i++) {
        for (int j = 0; j < stride; j++) {
            double count = (double) dye_tracks[i]->value->counts[j];
            raw_dataset[i * stride + j] = count;
        }
    }
    dataset = new Matrix<double>(raw_dataset, num_train, stride);
    index = new Index<L2<double>>(*dataset,
                                  AutotunedIndexParams(
                                          0.9f,  // target precision
                                          1.0f,  // build weight
                                          0.0f,  // memory weight
                                          0.0001f));  // sample fraction
    index->buildIndex();
}

KWANNClassifier::~KWANNClassifier() {
    delete[] dataset->ptr();
    delete dataset;
    delete index;
}

ScoredClassification KWANNClassifier::classify(const Radiometry& radiometry) {
    Matrix<double> query(radiometry.intensities,
                         1,  // num rows (num queries)
                         num_timesteps * num_channels);  // num columns
    Matrix<int> indices(new int[k],
                        1,  // num rows (num queries)
                        k);  // num results per query
    Matrix<double> dists_sq(new double[k],
                            1,  // num rows (num queries)
                            k);  // num results per query
    index->knnSearch(query,
                     indices,
                     dists_sq,
                     k,
                     SearchParams(FLANN_CHECKS_AUTOTUNED));
    delete[] dists_sq.ptr();
    unordered_map<int, double> id_score_map;
    for (int i = 0; i < k; i++) {
        int index = indices[0][i];
        SourcedData<DyeTrack*,
                    SourceCountMap<int>*>* dye_track = dye_tracks[index];
        double weight = 1.0;
        for (int j = 0; j < num_timesteps * num_channels; j++) {
            weight *= pdf(radiometry.intensities[j],
                          dye_track->value->counts[j]);
        }
        for (int j = 0; j < dye_track->source->num_sources; j++) {
            int id = dye_track->source->sources_with_counts[j]->source;
            int count = dye_track->source->sources_with_counts[j]->count;
            id_score_map[id] += weight * count;
        }
    }
    delete[] indices.ptr();
    int best_id = -1;
    double best_score = -1.0;
    double total_score = 0.0;
    for (const auto& id_and_score : id_score_map) {
        double id = id_and_score.first;
        double score = id_and_score.second;
        total_score += score;
        if (score > best_score) {
            best_id = id;
            best_score = score;
        }
    }
    return ScoredClassification(best_id, best_score, total_score);
}

ScoredClassification* KWANNClassifier::classify(int num_radiometries, 
                                                Radiometry** radiometries) {
    ScoredClassification* results = new ScoredClassification[num_radiometries];
    for (int i = 0; i < num_radiometries; i++) {
        results[i] = classify(*radiometries[i]);
    }
    return results;
}

}  // namespace fluoroseq
