/******************************************************************************\
* Author: Matthew Beauregard Smith                                             *
* Affiliation: The University of Texas at Austin                               *
* Department: Oden Institute and Institute for Cellular and Molecular Biology  *
* PI: Edward Marcotte                                                          *
* Project: Protein Fluorosequencing                                            *
\******************************************************************************/

#ifndef KD_TREE_K_BEST_H
#define KD_TREE_K_BEST_H

// Standard C++ library headers:
#include <cfloat>
#include <queue>
#include <utility>
#include <vector>

namespace fluoroseq {
namespace kd_tree {

template <typename T>
class KBest {
  public:
    KBest(int k) : k(k), kth_distance(DBL_MAX) {}

    virtual void consider(double d, T* t) {
        if (d < kth_distance) {
            if (pq.size() == k) {
                pq.pop();
            }
            pq.push(std::pair<double, T*>(d, t));
            if (pq.size() == k) {
                kth_distance = pq.top().first;
            }
        }
    }

    virtual void fill(std::vector<T*>* k_nearest,
                      std::vector<double>* dists_sq) {
        k_nearest->reserve(pq.size());
        dists_sq->reserve(pq.size());
        while (!pq.empty()) {
            dists_sq->push_back(pq.top().first);
            k_nearest->push_back(pq.top().second);
            pq.pop();
        }
    }

    int k;
    double kth_distance;
    std::priority_queue<std::pair<double, T*>> pq;
};

}  // namespace kd_tree
}  // namespace fluoroseq

#endif  // KD_TREE_K_BEST_H
