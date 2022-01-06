/******************************************************************************\
* Author: Matthew Beauregard Smith                                             *
* Affiliation: The University of Texas at Austin                               *
* Department: Oden Institute and Institute for Cellular and Molecular Biology  *
* PI: Edward Marcotte                                                          *
* Project: Protein Fluorosequencing                                            *
\******************************************************************************/

#ifndef WHATPROT_TENSOR_TENSOR_VECTOR_ITERATOR_H
#define WHATPROT_TENSOR_TENSOR_VECTOR_ITERATOR_H

// Local project headers:
#include "tensor/base-tensor-vector-iterator.h"

namespace whatprot {

class TensorVectorIterator : public BaseTensorVectorIterator<false> {
public:
    // This next line inherits all constructors of base class.
    using BaseTensorVectorIterator<false>::BaseTensorVectorIterator;
};

}  // namespace whatprot

#endif  // WHATPROT_TENSOR_TENSOR_VECTOR_ITERATOR_H
