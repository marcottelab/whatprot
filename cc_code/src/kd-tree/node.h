/******************************************************************************\
* Author: Matthew Beauregard Smith                                             *
* Affiliation: The University of Texas at Austin                               *
* Department: Oden Institute and Institute for Cellular and Molecular Biology  *
* PI: Edward Marcotte                                                          *
* Project: Protein Fluorosequencing                                            *
\******************************************************************************/

#ifndef KD_TREE_NODE_H
#define KD_TREE_NODE_H

// Standard C++ library headers:
#include <vector>

// Local project headers:
#include "kd-tree/k-best.h"

namespace fluoroseq {
namespace kd_tree {

template <typename E, typename Q>
class Node {
  public:
    virtual void search(const Q& query, KBest<E>* k_best) const = 0;
    virtual ~Node(){};
};

}  // namespace kd_tree
}  // namespace fluoroseq

#endif  // KD_TREE_NODE_H
