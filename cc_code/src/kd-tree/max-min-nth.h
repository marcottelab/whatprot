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

template <typename T>
T* partition(T* begin, T* end, int s, T* pivot_ptr, int* num_eq_to_pivot) {
    *num_eq_to_pivot = 0;
    T* left_ptr = begin;
    T* right_ptr = end - 1;
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

template <typename T>
void partition_alternate(T* begin, T* end, int s) {
    T* left_ptr = begin;
    T* right_ptr = end - 1;
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
template <typename T>
T* nth_element(T* begin, T* end, int s, T* nth) {
    int num_eq_to_pivot;
    T* pivot_ptr = partition(begin, end, s, nth, &num_eq_to_pivot);
    T* right_ptr = pivot_ptr;
    T* left_ptr = pivot_ptr - num_eq_to_pivot;
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

template <typename T>
double max_element(T* begin, T* end, int s) {
    double max = DBL_MIN;
    for (T* ptr = begin; ptr < end; ptr++) {
        max = std::max(max, (*ptr)[s]);
    }
    return max;
}

template <typename T>
double min_element(T* begin, T* end, int s) {
    double min = DBL_MAX;
    for (T* ptr = begin; ptr < end; ptr++) {
        min = std::min(min, (*ptr)[s]);
    }
    return min;
}

}  // namespace kd_tree
}  // namespace fluoroseq

#endif  // KD_TREE_MAX_MIN_NTH_H
