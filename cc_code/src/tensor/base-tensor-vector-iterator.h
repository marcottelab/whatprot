/******************************************************************************\
* Author: Matthew Beauregard Smith                                             *
* Affiliation: The University of Texas at Austin                               *
* Department: Oden Institute and Institute for Cellular and Molecular Biology  *
* PI: Edward Marcotte                                                          *
* Project: Protein Fluorosequencing                                            *
\******************************************************************************/

#ifndef WHATPROT_TENSOR_BASE_TENSOR_VECTOR_ITERATOR_H
#define WHATPROT_TENSOR_BASE_TENSOR_VECTOR_ITERATOR_H

// Standard C++ library headers:
#include <type_traits>

// Local project headers:
#include "tensor/tensor-iterator.h"
#include "tensor/vector.h"
#include "util/kd-range.h"

namespace whatprot {

// Templatizing constness to share code between const and non-const iterators.
template <bool is_const>
class BaseTensorVectorIterator {
public:
    BaseTensorVectorIterator(unsigned int order,
                             const KDRange& itr_range,
                             const KDRange& tsr_range,
                             const int* stride,
                             unsigned int size,
                             double* values,
                             unsigned int vector_dimension)
            : vector_min(tsr_range.min[vector_dimension]),
              vector_max(tsr_range.max[vector_dimension]),
              vector_stride(stride[vector_dimension]),
              modified_range(itr_range),
              iter(order, modified_range, tsr_range, size, values) {
        modified_range.min[vector_dimension] = tsr_range.min[vector_dimension];
        modified_range.max[vector_dimension] =
                1 + tsr_range.min[vector_dimension];
        iter.reset();
    }

    void reset() {
        iter.reset();
    }

    void advance() {
        iter.advance();
    }

    typename std::conditional<is_const, const Vector*, Vector*>::type get() {
        return new Vector(vector_min, vector_max, vector_stride, iter.get());
    }

    bool done() {
        return iter.done();
    }

    unsigned int vector_min;
    unsigned int vector_max;
    int vector_stride;
    KDRange modified_range;

    // iter is not a ConstTensorIterator, even if is_const is true. This may
    // seem a little strange, but we trust users of this class to not directly
    // access iter, even though it is public. A non-const tensor iterator is
    // necessary to avoid the need for a ConstVector type to implement the get()
    // function for the (is_const == true) case.
    TensorIterator iter;
};

}  // namespace whatprot

#endif  // WHATPROT_TENSOR_BASE_TENSOR_VECTOR_ITERATOR_H