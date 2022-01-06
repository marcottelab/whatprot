/******************************************************************************\
* Author: Matthew Beauregard Smith                                             *
* Affiliation: The University of Texas at Austin                               *
* Department: Oden Institute and Institute for Cellular and Molecular Biology  *
* PI: Edward Marcotte                                                          *
* Project: Protein Fluorosequencing                                            *
\******************************************************************************/

#ifndef WHATPROT_TENSOR_CONST_TENSOR_ITERATOR_H
#define WHATPROT_TENSOR_CONST_TENSOR_ITERATOR_H

// Local project headers:
#include "tensor/base-tensor-iterator.h"

namespace whatprot {

class ConstTensorIterator : public BaseTensorIterator<true> {
public:
    // This next line inherits all constructors of base class.
    using BaseTensorIterator<true>::BaseTensorIterator;
};

}  // namespace whatprot

#endif  // WHATPROT_TENSOR_CONST_TENSOR_ITERATOR_H
