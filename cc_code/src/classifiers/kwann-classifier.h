/******************************************************************************\
* Author: Matthew Beauregard Smith                                             *
* Affiliation: The University of Texas at Austin                               *
* Department: Oden Institute and Institute for Cellular and Molecular Biology  *
* PI: Edward Marcotte                                                          *
* Project: Protein Fluorosequencing                                            *
\******************************************************************************/

#ifndef FLUOROSEQ_CLASSIFIERS_KWANN_CLASSIFIER_H
#define FLUOROSEQ_CLASSIFIERS_KWANN_CLASSIFIER_H

// Standard C++ library headers:
#include <functional>
#include <unordered_map>
#include <vector>

// Local project headers:
#include "common/dye-track.h"
#include "common/error-model.h"
#include "common/radiometry.h"
#include "common/scored-classification.h"
#include "common/sourced-data.h"
#include "flann/flann.hpp"

namespace fluoroseq {

// KWANN stands for K Weighted, Approximate, Nearest-Neighbors. The weights come
// from the pdf (Probability Density Function) provided to the constructor. To
// find approximate nearest neighbors we use the FLANN (Fast Linear Approximate
// Nearest Neighbors) library (see extern/flann-*).
class KWANNClassifier {
  public:
    KWANNClassifier(
            int num_timesteps,
            int num_channels,
            const ErrorModel& error_model,
            int k,
            double sigma,
            const std::vector<SourcedData<DyeTrack, SourceCountHitsList<int>>>&
                    dye_tracks);
    ~KWANNClassifier();
    double classify_helper(const Radiometry& radiometry,
                           std::unordered_map<int, double>* id_score_map);
    ScoredClassification classify(const Radiometry& radiometry);
    std::vector<ScoredClassification> classify(const Radiometry& radiometry,
                                               int h);
    std::vector<ScoredClassification> classify(
            const std::vector<Radiometry>& radiometries);

    std::function<double(double, int)> kernel;
    flann::Index<flann::L2<double>> index;
    flann::Matrix<double> dataset;
    const std::vector<SourcedData<DyeTrack, SourceCountHitsList<int>>>&
            dye_tracks;
    int num_train;
    int num_timesteps;
    int num_channels;
    int k;  // number of approximate nearest neighbors to use
};

}  // namespace fluoroseq

#endif  // FLUOROSEQ_CLASSIFIERS_KWANN_CLASSIFIER_H