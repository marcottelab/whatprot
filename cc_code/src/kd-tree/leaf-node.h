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
            consider(query, t, k_best);
        }
    }

    void consider(const Q& query, E* entry, KBest<E>* k_best) const {
        double kth_dist_sq = k_best->kth_dist_sq;
        double dist_sq = 0.0;
        // This is called loop unrolling. We stretch out a loop in order to
        // check the loop condition less frequently, which is more efficient.
        // This is particularly important because of the (dist_sq > kth_dist_sq)
        // condition, which is expensive because the computation can't start
        // until the computation for dist_sq is finished.
        int i = 0;
        while (i < d - 3) {
            double x1 = query[i] - (*entry)[i];
            double x2 = query[i + 1] - (*entry)[i + 1];
            double x3 = query[i + 2] - (*entry)[i + 2];
            double x4 = query[i + 3] - (*entry)[i + 3];
            dist_sq += x1 * x1 + x2 * x2 + x3 * x3 + x4 * x4;
            // Here we consider returning early to avoid unnecessary additional
            // calculation.
            if (dist_sq >= kth_dist_sq) {
                return;
            }
            i += 4;
        }
        while (i < d) {
            double x = query[i] - (*entry)[i];
            dist_sq += x * x;
            i++;
        }
        if (dist_sq >= kth_dist_sq) {
            return;
        }
        k_best->insert(dist_sq, entry);
    }

    int d;
    E* begin;
    E* end;
};

}  // namespace kd_tree
}  // namespace fluoroseq

#endif  // KD_TREE_LEAF_NODE_H
