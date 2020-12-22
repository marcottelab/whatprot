/******************************************************************************\
* Author: Matthew Beauregard Smith                                             *
* Affiliation: The University of Texas at Austin                               *
* Department: Oden Institute and Institute for Cellular and Molecular Biology  *
* PI: Edward Marcotte                                                          *
* Project: Protein Fluorosequencing                                            *
\******************************************************************************/

#ifndef FLUOROSEQ_TENSOR_TENSOR_H
#define FLUOROSEQ_TENSOR_TENSOR_H

// Local project headers:
#include "tensor/tensor-iterator.h"

namespace fluoroseq {

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

}  // namespace fluoroseq

#endif  // FLUOROSEQ_TENSOR_TENSOR_H
