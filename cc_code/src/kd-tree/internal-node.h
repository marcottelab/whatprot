/******************************************************************************\
* Author: Matthew Beauregard Smith                                             *
* Affiliation: The University of Texas at Austin                               *
* Department: Oden Institute and Institute for Cellular and Molecular Biology  *
* PI: Edward Marcotte                                                          *
* Project: Protein Fluorosequencing                                            *
\******************************************************************************/

#ifndef KD_TREE_INTERNAL_NODE_H
#define KD_TREE_INTERNAL_NODE_H

// Local project headers:
#include "kd-tree/k-best.h"
#include "kd-tree/node.h"

namespace fluoroseq {
namespace kd_tree {

template <typename E, typename Q>
class InternalNode : public Node<E, Q> {
public:
    InternalNode(Node<E, Q>* left_child,
                 Node<E, Q>* right_child,
                 double max_left,
                 double min_right,
                 int s)
            : left_child(left_child),
              right_child(right_child),
              max_left(max_left),
              min_right(min_right),
              split_value((max_left + min_right) / 2),
              s(s) {}

    virtual ~InternalNode() {
        delete left_child;
        delete right_child;
    }

    virtual void search(const Q& query, KBest<E>* k_best) const {
        double query_value = query[s];
        if (query_value < split_value) {
            left_child->search(query, k_best);
            // You CANNOT just compare the dist directly for the second search.
            // These are squared distances, for performance reasons. BE CAREFUL
            // IF EDITING THIS CODE.
            double right_dist = min_right - query_value;
            double right_dist_sq = right_dist * right_dist;
            if (k_best->kth_dist_sq > right_dist_sq) {
                right_child->search(query, k_best);
            }
        } else {
            right_child->search(query, k_best);
            // Again, must use squared distances. See previous comment.
            double left_dist = query_value - max_left;
            double left_dist_sq = left_dist * left_dist;
            if (k_best->kth_dist_sq > left_dist_sq) {
                left_child->search(query, k_best);
            }
        }
    }

    Node<E, Q>* left_child;
    Node<E, Q>* right_child;
    double max_left;
    double min_right;
    double split_value;
    int s;
};

}  // namespace kd_tree
}  // namespace fluoroseq

#endif  // KD_TREE_INTERNAL_NODE_H
