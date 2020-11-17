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
    // IMPORTANT: There is ONE TensorIterator. Never try to hold multiple
    // TensorIterators on the same Tensor at the same time. Furthermore, the
    // caller of the iterator() function DOES NOT OWN the TensorIterator that is
    // returned.
    TensorIterator* iterator();

    TensorIterator* tensor_iterator;
    double* values;
    int* shape;
    int* strides;
    int size;
    int order;
};

}  // namespace fluoroseq

#endif  // FLUOROSEQ_TENSOR_TENSOR_H
