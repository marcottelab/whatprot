/******************************************************************************\
* Author: Matthew Beauregard Smith                                             *
* Affiliation: The University of Texas at Austin                               *
* Department: Oden Institute and Institute for Cellular and Molecular Biology  *
* PI: Edward Marcotte                                                          *
* Project: Protein Fluorosequencing                                            *
\******************************************************************************/

#ifndef WHATPROT_TENSOR_TENSOR_H
#define WHATPROT_TENSOR_TENSOR_H

// Local project headers:
#include "tensor/tensor-iterator.h"

namespace whatprot {

class Tensor {
public:
    Tensor(int order, int* shape);
    Tensor(Tensor&& other);
    ~Tensor();
    double& operator[](int* loc);
    TensorIterator* iterator();

    double* values;
    int* shape;
    int* strides;
    int size;
    int order;
};

}  // namespace whatprot

#endif  // WHATPROT_TENSOR_TENSOR_H
