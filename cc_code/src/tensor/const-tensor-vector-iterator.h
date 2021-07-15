/******************************************************************************\
* Author: Matthew Beauregard Smith                                             *
* Affiliation: The University of Texas at Austin                               *
* Department: Oden Institute and Institute for Cellular and Molecular Biology  *
* PI: Edward Marcotte                                                          *
* Project: Protein Fluorosequencing                                            *
\******************************************************************************/

#ifndef WHATPROT_TENSOR_CONST_TENSOR_VECTOR_ITERATOR_H
#define WHATPROT_TENSOR_CONST_TENSOR_VECTOR_ITERATOR_H

// Local project headers:
#include "tensor/base-tensor-vector-iterator.h"

namespace whatprot {

class ConstTensorVectorIterator : public BaseTensorVectorIterator<true> {
public:
    // This next line inherits all constructors of base class.
    using BaseTensorVectorIterator<true>::BaseTensorVectorIterator;
};

}  // namespace whatprot

#endif  // WHATPROT_TENSOR_CONST_TENSOR_VECTOR_ITERATOR_H
