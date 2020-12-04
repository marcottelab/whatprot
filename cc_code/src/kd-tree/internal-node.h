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

template <typename T>
class InternalNode : public Node<T> {
  public:
    InternalNode(Node<T>* left_child,
                 Node<T>* right_child,
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

    virtual void search(const T& query, KBest<T>* k_best) const {
        double query_value = query[s];
        if (query_value < split_value) {
            left_child->search(query, k_best);
            if (query_value + k_best->kth_distance > min_right) {
                right_child->search(query, k_best);
            }
        } else {
            right_child->search(query, k_best);
            if (query_value - k_best->kth_distance < max_left) {
                left_child->search(query, k_best);
            }
        }
    }

    Node<T>* left_child;
    Node<T>* right_child;
    double max_left;
    double min_right;
    double split_value;
    int s;
};

}  // namespace kd_tree
}  // namespace fluoroseq

#endif  // KD_TREE_INTERNAL_NODE_H
