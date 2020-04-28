// Author: Matthew Beauregard Smith (UT Austin)
#include "tensor.h"

#include <algorithm>  // needed for std::copy

#include "tensor/tensor_iterator.h"

namespace fluoroseq {

namespace {
using std::copy;
}

Tensor::Tensor(int order, int* shape) : order(order) {
    this->shape = new int[order];
    copy(shape, shape + order, this->shape);
    size = 1;
    strides = new int[order];
    for (int i = order - 1; i >= 0; i--) {
        strides[i] = size;
        size *= shape[i];
    }
    values = new double[size];
}

Tensor::~Tensor() {
    delete[] shape;
    delete[] strides;
    delete[] values;
}

double& Tensor::operator[](int* loc) {
    int index = 0;
    for (int i = 0; i < order; i++) {
        index += strides[i] * loc[i];
    }
    return values[index];
}

TensorIterator* Tensor::iterator() {
    return new TensorIterator(order, shape, size, values);
}

}  // namespace fluoroseq
