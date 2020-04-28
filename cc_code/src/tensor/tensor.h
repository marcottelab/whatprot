// Author: Matthew Beauregard Smith (UT Austin)
#ifndef FLUOROSEQ_TENSOR_TENSOR_H
#define FLUOROSEQ_TENSOR_TENSOR_H

#include "tensor/tensor_iterator.h"

namespace fluoroseq {

class Tensor {
public:
    Tensor(int order, int* shape);
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
