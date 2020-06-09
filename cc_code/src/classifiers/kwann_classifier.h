// Author: Matthew Beauregard Smith (UT Austin)
#ifndef FLUOROSEQ_CLASSIFIERS_KWANN_CLASSIFIER_H
#define FLUOROSEQ_CLASSIFIERS_KWANN_CLASSIFIER_H

#include <functional>
#include <unordered_map>

#include "common/dye_track.h"
#include "common/radiometry.h"
#include "common/scored_classification.h"
#include "common/sourced_data.h"
#include "flann/flann.hpp"

namespace fluoroseq {

// KWANN stands for K Weighted, Approximate, Nearest-Neighbors. The weights come
// from the pdf (Probability Density Function) provided to the constructor. To
// find approximate nearest neighbors we use the FLANN (Fast Linear Approximate
// Nearest Neighbors) library (see extern/flann-*).
class KWANNClassifier {
public:
    KWANNClassifier(int num_timesteps,
                    int num_channels,
                    std::function<double (double, int)> pdf,
                    int k,
                    int num_train,
                    SourcedData<DyeTrack*, SourceCountMap<int>*>** dye_tracks);
    ~KWANNClassifier();
    ScoredClassification classify(const Radiometry& radiometry);
    ScoredClassification* classify(int num_radiometries, 
                                    Radiometry** radiometries);

    std::function<double (double, int)> pdf;
    flann::Index<flann::L2<double>>* index;
    flann::Matrix<double>* dataset;
    SourcedData<DyeTrack*, SourceCountMap<int>*>** dye_tracks;  // not owned
    int num_train;
    int num_timesteps;
    int num_channels;
    int k;  // number of approximate nearest neighbors to use
};

}  // namespace fluoroseq

#endif  // FLUOROSEQ_CLASSIFIERS_KWANN_CLASSIFIER_H