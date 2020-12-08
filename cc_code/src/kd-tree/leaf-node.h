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

template <typename E, typename Q>
class LeafNode : public Node<E, Q> {
  public:
    LeafNode(int d, E* begin, E* end) : d(d), begin(begin), end(end) {}
    virtual ~LeafNode() {}

    virtual void search(const Q& query, KBest<E>* k_best) const {
        for (E* t = begin; t < end; t++) {
            double dist = distance(query, *t);
            k_best->consider(dist, t);
        }
    }

    double distance(const Q& query, const E& entry) const {
        double dist = 0.0;
        for (int i = 0; i < d; i++) {
            double diff = query[i] - entry[i];
            dist += diff * diff;
        }
        return dist;
    }

    int d;
    E* begin;
    E* end;
};

}  // namespace kd_tree
}  // namespace fluoroseq

#endif  // KD_TREE_LEAF_NODE_H
