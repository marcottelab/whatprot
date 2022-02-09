/******************************************************************************\
* Author: Matthew Beauregard Smith                                             *
* Affiliation: The University of Texas at Austin                               *
* Department: Oden Institute and Institute for Cellular and Molecular Biology  *
* PI: Edward Marcotte                                                          *
* Project: Protein Fluorosequencing                                            *
\******************************************************************************/

#ifndef WHATPROT_TENSOR_BASE_TENSOR_ITERATOR_H
#define WHATPROT_TENSOR_BASE_TENSOR_ITERATOR_H

// Standard C++ library headers:
#include <type_traits>

// Local project headers:
#include "util/kd-range.h"

namespace whatprot {

// Templatizing constness to share code between const and non-const iterators.
template <bool is_const>
class BaseTensorIterator {
public:
    BaseTensorIterator(
            unsigned int order,
            const KDRange& itr_range,
            const KDRange& tsr_range,
            unsigned int size,
            typename std::conditional<is_const, const double*, double*>::type
                    values)
            : values(values),
              itr_range(itr_range),
              tsr_range(tsr_range),
              order(order),
              index(0),
              size(size),
              is_done(itr_range.is_empty()) {
        loc = new unsigned int[order]();
        reset();
    }

    ~BaseTensorIterator() {
        delete[] loc;
    }

    void reset() {
        index = 0;
        unsigned int stride = 1;  // step size for current value of o.
        // Need to use signed int o to detect when out of entries to decrease.
        for (int o = order - 1; o >= 0; o--) {
            loc[o] = itr_range.min[o];
            index += stride * (itr_range.min[o] - tsr_range.min[o]);
            stride *= tsr_range.max[o] - tsr_range.min[o];
        }
    }

    void set_to_last() {
        index = 0;
        unsigned int stride = 1;  // step size for current value of o.
        // Need to use signed int o to detect when out of entries to decrease.
        for (int o = order - 1; o >= 0; o--) {
            loc[o] = itr_range.max[o] - 1;
            index += stride * (itr_range.max[o] - tsr_range.min[o] - 1);
            stride *= tsr_range.max[o] - tsr_range.min[o];
        }
    }

    void advance() {
        index++;
        unsigned int stride = 1;  // step size for current value of o.
        // Need to use signed int o to detect when out of entries to decrease.
        for (int o = order - 1; o >= 0; o--) {
            loc[o]++;
            if (loc[o] < itr_range.max[o]) {
                return;
            } else {
                // (itr_range.min[o] - tsr_range.min[o] + tsr_range.max[o]
                // - itr_range.max[o]) is the number of entries for this order
                // between the max of one set of numbers in range to the min of
                // the next.
                index += stride
                         * (itr_range.min[o] - tsr_range.min[o]
                            + tsr_range.max[o] - itr_range.max[o]);
                // Stride needs to be raised to account for one more dimension.
                stride *= tsr_range.max[o] - tsr_range.min[o];
                loc[o] = itr_range.min[o];
            }
        }
        // Can only get here if we reset every order of the tensor, which means
        // we're done.
        is_done = true;
    }

    void retreat() {
        index--;
        unsigned int stride = 1;  // step size for current value of o.
        // Need to use signed int o to detect when out of entries to decrease.
        for (int o = order - 1; o >= 0; o--) {
            if (loc[o] > itr_range.min[o]) {
                loc[o]--;
                return;
            } else {
                // (itr_range.min[o] - tsr_range.min[o] + tsr_range.max[o]
                // - itr_range.max[o]) is the number of entries for this order
                // between the max of one set of numbers in range to the min of
                // the next.
                index -= stride
                         * (itr_range.min[o] - tsr_range.min[o]
                            + tsr_range.max[o] - itr_range.max[o]);
                // Stride needs to be raised to account for one more dimension.
                stride *= tsr_range.max[o] - tsr_range.min[o];
                loc[o] = itr_range.max[o] - 1;
            }
        }
        // Can only get here if we reset every order of the tensor, which means
        // we're done.
        is_done = true;
    }

    typename std::conditional<is_const, const double*, double*>::type get() {
        return &values[index];
    }

    bool done() {
        return is_done;
    }

    typename std::conditional<is_const, const double*, double*>::type
            values;  // not owned
    const KDRange& itr_range;  // not owned
    const KDRange& tsr_range;  // not owned
    unsigned int* loc;
    const unsigned int order;
    unsigned int index;  // current index directly into values
    const unsigned int size;  // length of values
    bool is_done;
};

}  // namespace whatprot

#endif  // WHATPROT_TENSOR_BASE_TENSOR_ITERATOR_H