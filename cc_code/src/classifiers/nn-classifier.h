/******************************************************************************\
* Author: Matthew Beauregard Smith                                             *
* Affiliation: The University of Texas at Austin                               *
* Department: Oden Institute and Institute for Cellular and Molecular Biology  *
* PI: Edward Marcotte                                                          *
* Project: Protein Fluorosequencing                                            *
\******************************************************************************/

#ifndef WHATPROT_CLASSIFIERS_NN_CLASSIFIER_H
#define WHATPROT_CLASSIFIERS_NN_CLASSIFIER_H

// Standard C++ library headers:
#include <functional>
#include <unordered_map>
#include <utility>  // so that we can overload the std::swap function.
#include <vector>

// Local project headers:
#include "common/dye-track.h"
#include "common/radiometry.h"
#include "common/scored-classification.h"
#include "common/sourced-data.h"
#include "kd-tree/kd-tree.h"
#include "parameterization/model/sequencing-model.h"

namespace whatprot {

class KDTEntry {
public:
    KDTEntry(SourcedData<DyeTrack, SourceCountHitsList<int>>&& dye_track);
    KDTEntry(KDTEntry&& other);
    KDTEntry& operator=(KDTEntry&& other);
    double operator[](int i) const;
    SourcedData<DyeTrack, SourceCountHitsList<int>> dye_track;
    int hits;
};

}  // namespace whatprot

// Namespace injection to make KDTEntry swappable.
namespace std {
void swap(whatprot::KDTEntry& e1, whatprot::KDTEntry& e2);
}  // namespace std

namespace whatprot {

class KDTQuery {
public:
    KDTQuery(const Radiometry& rad);
    double operator[](int i) const;
    Radiometry rad;
};

class NNClassifier {
public:
    NNClassifier(unsigned int num_timesteps,
                 unsigned int num_channels,
                 const SequencingModel& seq_model,
                 int k,
                 double sig,
                 std::vector<SourcedData<DyeTrack, SourceCountHitsList<int>>>*
                         dye_tracks);
    ~NNClassifier();
    double classify_helper(const Radiometry& radiometry,
                           std::unordered_map<int, double>* id_score_map);
    ScoredClassification classify(const Radiometry& radiometry);
    std::vector<ScoredClassification> classify(const Radiometry& radiometry,
                                               unsigned int h);
    std::vector<ScoredClassification> classify(
            const std::vector<Radiometry>& radiometries);

    KDTree<KDTEntry, KDTQuery>* kd_tree;
    int num_train;
    unsigned int num_timesteps;
    unsigned int num_channels;
    int k;  // number of nearest neighbors to use
    double two_sig_sq;  // sig to use for kernel weighting
};

}  // namespace whatprot

#endif  // WHATPROT_CLASSIFIERS_NN_CLASSIFIER_H