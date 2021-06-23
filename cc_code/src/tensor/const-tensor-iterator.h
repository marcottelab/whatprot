/******************************************************************************\
* Author: Matthew Beauregard Smith                                             *
* Affiliation: The University of Texas at Austin                               *
* Department: Oden Institute and Institute for Cellular and Molecular Biology  *
* PI: Edward Marcotte                                                          *
* Project: Protein Fluorosequencing                                            *
\******************************************************************************/

#ifndef WHATPROT_TENSOR_CONST_TENSOR_ITERATOR_H
#define WHATPROT_TENSOR_CONST_TENSOR_ITERATOR_H

#include "tensor/base-tensor-iterator.h"

namespace whatprot {

class ConstTensorIterator : public BaseTensorIterator<const double*> {
public:
    // This next line inherits all constructors of base class.
    using BaseTensorIterator<const double*>::BaseTensorIterator;
};

}  // namespace whatprot

#endif  // WHATPROT_TENSOR_CONST_TENSOR_ITERATOR_H
