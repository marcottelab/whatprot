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

// Template requirements for KDTree to work.
//   * typename E must have an operator[] function which receives an int and
//     returns a double.
//   * typename E must have a int member called size indicating the number of
//     simulated elements it represents.
//   * typename E must be able to be swapped. std::swap(E, E) needs to work.
//     For the best performance, consider doing this through move operations.
//   * typename Q must have an operator[] function which receives an int and
//     returns a double.
template <typename E, typename Q>
class KDTree {
public:
    // At construction, you must give the following.
    //   1. k - the number of nearest neighbors to be returned when a search is
    //      performed.
    //   2. d - the dimensionality of your data.
    //   3. v - a vector of E typed elements, which must give a valid response
    //      when called with the operator[] function for any element in the
    //      range [0, d).
    KDTree(int k, int d, std::vector<E>&& v)
            : k(k),
              values(std::move(v)),
              root(kd_tree::make_node<E, Q>(
                      k, d, &values[0], &values[values.size()])) {}

    ~KDTree() {
        delete root;
    }

    void search(const Q& query,
                std::vector<E*>* k_nearest,
                std::vector<double>* dists_sq) const {
        kd_tree::KBest<E> k_best(k);
        root->search(query, &k_best);
        k_best.fill(k_nearest, dists_sq);
    }

    int k;
    std::vector<E> values;
    kd_tree::Node<E, Q>* root;
};

}  // namespace fluoroseq

#endif  // KD_TREE_KD_TREE_H
