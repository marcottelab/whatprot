/******************************************************************************\
* Author: Matthew Beauregard Smith                                             *
* Affiliation: The University of Texas at Austin                               *
* Department: Oden Institute and Institute for Cellular and Molecular Biology  *
* PI: Edward Marcotte                                                          *
* Project: Protein Fluorosequencing                                            *
\******************************************************************************/

#ifndef KD_TREE_MAKE_NODE_H
#define KD_TREE_MAKE_NODE_H

// Standard C++ library headers:
#include <algorithm>
#include <cfloat>
#include <functional>
#include <vector>

// Local project headers:
#include "kd-tree/internal-node.h"
#include "kd-tree/leaf-node.h"
#include "kd-tree/max-min-nth.h"
#include "kd-tree/node.h"

namespace fluoroseq {
namespace kd_tree {

template <typename E, typename Q>
Node<E, Q>* make_node(int k, int d, E* begin, E* end) {
    // We want to split on the dimension with the biggest range. To do that, we
    // find the minimums and maximums of every dimension, which we can later use
    // to find the range of every dimension and take the max.
    std::vector<double> mins(d, DBL_MAX);
    std::vector<double> maxes(d, DBL_MIN);
    for (E* t = begin; t < end; t++) {
        for (int i = 0; i < d; i++) {
            mins[i] = std::min(mins[i], (*t)[i]);
            maxes[i] = std::max(maxes[i], (*t)[i]);
        }
    }
    // Take largest range, computing using mins and maxes.
    double best_range = -1.0;
    int s = -1;  // s is the split_dim
    for (int i = 0; i < d; i++) {
        double range = maxes[i] - mins[i];
        if (range > best_range) {
            best_range = range;
            s = i;
        }
    }
    // Now we find the middle using reference math with begin and end. This is
    // what we will use as the "nth_element" below.
    E* nth = &begin[(size_t)(end - begin) / 2];
    // nth_element()
    //   * guarantees that every element to the left of the nth element will be
    //     strictly less the nth element and every element to its right.
    //   * However, the exact location of the nth element may change - the
    //     algorithm will do its best to change this location as little as
    //     possible (this is necessary to deal with duplicate values).
    nth = nth_element(begin, end, s, nth);
    // Set ranges for left and right children and create them. Important
    // points:
    //   * "end" always signifies one PAST the last element, therefore the end
    //     for the left node is the same as begin for the right node.
    //   * We don't want any children smaller than k. We plan to search for k
    //     nearest neighbors, therefore we waste time with reference hops to
    //     deal with any leaf that has less than k values. If either child of
    //     the new node would be too small, we can just make a leaf instead.
    int left_size = (int)(nth - begin);
    int right_size = (int)(end - nth);
    if (left_size < k || right_size < k) {
        return new LeafNode<E, Q>(d, begin, end);
    } else {
        double max_left = max_element(begin, nth, s);
        double min_right = min_element(nth, end, s);
        Node<E, Q>* left_child = make_node<E, Q>(k, d, begin, nth);
        Node<E, Q>* right_child = make_node<E, Q>(k, d, nth, end);
        return new InternalNode<E, Q>(
                left_child, right_child, max_left, min_right, s);
    }
}

}  // namespace kd_tree
}  // namespace fluoroseq

#endif  // KD_TREE_MAKE_NODE_H
