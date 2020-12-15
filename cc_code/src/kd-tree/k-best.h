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

template <typename E>
class KBest {
public:
    KBest(int k) : k(k), hits(0), kth_dist_sq(DBL_MAX) {}

    virtual void insert(double d, E* entry) {
        hits += entry->hits;
        pq.push(std::pair<double, E*>(d, entry));
        // Here we pop the top element of pq if removing it still allows us to
        // have a hits larger than k.
        int top_hits = pq.top().second->hits;
        if (hits - top_hits >= k) {
            hits -= top_hits;
            pq.pop();
        }
        // We only want to reset kth_dist_sq if we have enough elements. This
        // ensures that the kth_dist_sq remains DBL_MAX so that everything tried
        // will be added until there are enough elements.
        if (hits >= k) {
            kth_dist_sq = pq.top().first;
        }
    }

    virtual void fill(std::vector<E*>* k_nearest,
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
    int hits;
    double kth_dist_sq;
    std::priority_queue<std::pair<double, E*>> pq;
};

}  // namespace kd_tree
}  // namespace fluoroseq

#endif  // KD_TREE_K_BEST_H
