/******************************************************************************\
* Author: Matthew Beauregard Smith                                             *
* Affiliation: The University of Texas at Austin                               *
* Department: Oden Institute and Institute for Cellular and Molecular Biology  *
* PI: Edward Marcotte                                                          *
* Project: Protein Fluorosequencing                                            *
\******************************************************************************/

#ifndef KD_TREE_MAX_MIN_NTH_H
#define KD_TREE_MAX_MIN_NTH_H

// Standard C++ library headers:
#include <algorithm>
#include <cfloat>

namespace fluoroseq {
namespace kd_tree {

template <typename E>
E* partition(E* begin, E* end, int s, E* pivot_ptr, int* num_eq_to_pivot) {
    *num_eq_to_pivot = 0;
    E* left_ptr = begin;
    E* right_ptr = end - 1;
    double pivot_val = (*pivot_ptr)[s];
    std::swap(*pivot_ptr, *right_ptr);
    while (left_ptr < right_ptr) {
        while (left_ptr < right_ptr) {
            double left_val = (*left_ptr)[s];
            if (left_val == pivot_val) {
                (*num_eq_to_pivot)++;
            }
            if (left_val > pivot_val) {
                std::swap(*left_ptr, *right_ptr);
                right_ptr--;
                break;
            }
            left_ptr++;
        }
        while (left_ptr < right_ptr) {
            double right_val = (*right_ptr)[s];
            if (right_val == pivot_val) {
                (*num_eq_to_pivot)++;
            }
            if (right_val <= pivot_val) {
                std::swap(*right_ptr, *left_ptr);
                left_ptr++;
                break;
            }
            right_ptr--;
        }
    }
    return left_ptr;
}

template <typename E>
void partition_alternate(E* begin, E* end, int s) {
    E* left_ptr = begin;
    E* right_ptr = end - 1;
    double pivot_val = (*right_ptr)[s];
    while (left_ptr < right_ptr) {
        while (left_ptr < right_ptr) {
            double left_val = (*left_ptr)[s];
            if (left_val >= pivot_val) {
                std::swap(*left_ptr, *right_ptr);
                right_ptr--;
                break;
            }
            left_ptr++;
        }
        while (left_ptr < right_ptr) {
            double right_val = (*right_ptr)[s];
            if (right_val < pivot_val) {
                std::swap(*right_ptr, *left_ptr);
                left_ptr++;
                break;
            }
            right_ptr--;
        }
    }
}

// nth_element() will reorder elements from begin to end. Once run, it
// guarantees that elements in [begin, nth) are all strictly less than elements
// in [nth, end). Since there may be duplicates, nth_element() may need to move
// the pointer for nth. It will try to do this in the way that moves nth the
// least it can from its original position, while still satisfying the
// conditions.
template <typename E>
E* nth_element(E* begin, E* end, int s, E* nth) {
    int num_eq_to_pivot;
    E* pivot_ptr = partition(begin, end, s, nth, &num_eq_to_pivot);
    E* right_ptr = pivot_ptr;
    E* left_ptr = pivot_ptr - num_eq_to_pivot;
    if (left_ptr <= nth && nth <= right_ptr) {
        int num_pivots_to_right = right_ptr - nth;
        int num_pivots_to_left = nth - left_ptr;
        if (num_pivots_to_left > num_pivots_to_right) {
            nth = right_ptr + 1;
        } else {
            partition_alternate(begin, right_ptr + 1, s);
            nth = left_ptr;
        }
    } else if (nth < left_ptr) {
        partition_alternate(begin, right_ptr + 1, s);
        nth = nth_element(begin, left_ptr, s, nth);
    } else /* (nth > right_ptr) */ {
        nth = nth_element(right_ptr + 1, end, s, nth);
    }
    return nth;
}

template <typename E>
double max_element(E* begin, E* end, int s) {
    double max = DBL_MIN;
    for (E* ptr = begin; ptr < end; ptr++) {
        max = std::max(max, (*ptr)[s]);
    }
    return max;
}

template <typename E>
double min_element(E* begin, E* end, int s) {
    double min = DBL_MAX;
    for (E* ptr = begin; ptr < end; ptr++) {
        min = std::min(min, (*ptr)[s]);
    }
    return min;
}

}  // namespace kd_tree
}  // namespace fluoroseq

#endif  // KD_TREE_MAX_MIN_NTH_H
