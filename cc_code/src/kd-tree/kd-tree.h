/******************************************************************************\
* Author: Matthew Beauregard Smith                                             *
* Affiliation: The University of Texas at Austin                               *
* Department: Oden Institute and Institute for Cellular and Molecular Biology  *
* PI: Edward Marcotte                                                          *
* Project: Protein Fluorosequencing                                            *
\******************************************************************************/

#ifndef KD_TREE_KD_TREE_H
#define KD_TREE_KD_TREE_H

// Standard C++ library headers:
#include <functional>
#include <utility>
#include <vector>

// Local project headers:
#include "kd-tree/k-best.h"
#include "kd-tree/make-node.h"
#include "kd-tree/node.h"

namespace fluoroseq {

template <typename T>
class KDTree {
  public:
    KDTree(int k, int d, std::vector<T>&& v)
            : k(k),
              values(std::move(v)),
              root(kd_tree::make_node(
                      k, d, &values[0], &values[values.size()])) {}

    ~KDTree() {
        delete root;
    }

    void search(const T& query,
                std::vector<T*>* k_nearest,
                std::vector<double>* dists_sq) const {
        kd_tree::KBest<T> k_best(k);
        root->search(query, &k_best);
        k_best.fill(k_nearest, dists_sq);
    }

    int k;
    std::vector<T> values;
    kd_tree::Node<T>* root;
};

}  // namespace fluoroseq

#endif  // KD_TREE_KD_TREE_H
