/******************************************************************************\
* Author: Matthew Beauregard Smith                                             *
* Affiliation: The University of Texas at Austin                               *
* Department: Oden Institute and Institute for Cellular and Molecular Biology  *
* PI: Edward Marcotte                                                          *
* Project: Protein Fluorosequencing                                            *
\******************************************************************************/

#ifndef FLUOROSEQ_CLASSIFIERS_NN_CLASSIFIER_H
#define FLUOROSEQ_CLASSIFIERS_NN_CLASSIFIER_H

// Standard C++ library headers:
#include <functional>
#include <unordered_map>
#include <utility>  // so that we can overload the std::swap function.
#include <vector>

// Local project headers:
#include "common/dye-track.h"
#include "common/error-model.h"
#include "common/radiometry.h"
#include "common/scored-classification.h"
#include "common/sourced-data.h"
#include "kd-tree/kd-tree.h"

namespace fluoroseq {

class KDTEntry {
public:
    KDTEntry(SourcedData<DyeTrack, SourceCountHitsList<int>>&& dye_track);
    KDTEntry(KDTEntry&& other);
    KDTEntry& operator=(KDTEntry&& other);
    double operator[](int i) const;
    SourcedData<DyeTrack, SourceCountHitsList<int>> dye_track;
    int size;
};

}  // namespace fluoroseq

// Namespace injection to make KDTEntry swappable.
namespace std {
void swap(fluoroseq::KDTEntry& e1, fluoroseq::KDTEntry& e2);
}  // namespace std

namespace fluoroseq {

class KDTQuery {
public:
    KDTQuery(const Radiometry& rad);
    double operator[](int i) const;
    Radiometry rad;
};

class NNClassifier {
public:
    NNClassifier(int num_timesteps,
                 int num_channels,
                 const ErrorModel& error_model,
                 int k,
                 double sigma,
                 std::vector<SourcedData<DyeTrack, SourceCountHitsList<int>>>*
                         dye_tracks);
    ~NNClassifier();
    double classify_helper(const Radiometry& radiometry,
                           std::unordered_map<int, double>* id_score_map);
    ScoredClassification classify(const Radiometry& radiometry);
    std::vector<ScoredClassification> classify(const Radiometry& radiometry,
                                               int h);
    std::vector<ScoredClassification> classify(
            const std::vector<Radiometry>& radiometries);

    KDTree<KDTEntry, KDTQuery>* kd_tree;
    int num_train;
    int num_timesteps;
    int num_channels;
    int k;  // number of nearest neighbors to use
    double two_sigma_sq;  // sigma to use for kernel weighting
};

}  // namespace fluoroseq

#endif  // FLUOROSEQ_CLASSIFIERS_NN_CLASSIFIER_H