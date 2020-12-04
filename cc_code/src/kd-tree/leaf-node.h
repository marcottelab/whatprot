/******************************************************************************\
* Author: Matthew Beauregard Smith                                             *
* Affiliation: The University of Texas at Austin                               *
* Department: Oden Institute and Institute for Cellular and Molecular Biology  *
* PI: Edward Marcotte                                                          *
* Project: Protein Fluorosequencing                                            *
\******************************************************************************/

#ifndef KD_TREE_LEAF_NODE_H
#define KD_TREE_LEAF_NODE_H

// Local project headers:
#include "kd-tree/k-best.h"
#include "kd-tree/node.h"

namespace fluoroseq {
namespace kd_tree {

template <typename T>
class LeafNode : public Node<T> {
  public:
    LeafNode(int d, T* begin, T* end) : d(d), begin(begin), end(end) {}
    virtual ~LeafNode() {}

    virtual void search(const T& query, KBest<T>* k_best) const {
        for (T* t = begin; t < end; t++) {
            double dist = distance(query, *t);
            k_best->consider(dist, t);
        }
    }

    double distance(const T& t1, const T& t2) const {
        double dist = 0.0;
        for (int i = 0; i < d; i++) {
            double diff = t2[i] - t1[i];
            dist += diff * diff;
        }
        return dist;
    }

    int d;
    T* begin;
    T* end;
};

}  // namespace kd_tree
}  // namespace fluoroseq

#endif  // KD_TREE_LEAF_NODE_H
