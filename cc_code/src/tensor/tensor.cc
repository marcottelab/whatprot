/******************************************************************************\
* Author: Matthew Beauregard Smith                                             *
* Affiliation: The University of Texas at Austin                               *
* Department: Oden Institute and Institute for Cellular and Molecular Biology  *
* PI: Edward Marcotte                                                          *
* Project: Protein Fluorosequencing                                            *
\******************************************************************************/

// File under test:
#include "tensor.h"

// Standard C++ library headers:
#include <algorithm>  // needed for std::copy

// Local project headers:
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
    tensor_iterator = new TensorIterator(order, this->shape, size, values);
}

Tensor::Tensor(Tensor&& other)
        : tensor_iterator(other.tensor_iterator),
          values(other.values),
          shape(other.shape),
          strides(other.strides),
          size(other.size),
          order(other.order) {
    other.tensor_iterator = NULL;
    other.values = NULL;
    other.shape = NULL;
    other.strides = NULL;
}

Tensor::~Tensor() {
    if (shape != NULL) {
        delete[] shape;
    }
    if (strides != NULL) {
        delete[] strides;
    }
    if (values != NULL) {
        delete[] values;
    }
    if (tensor_iterator != NULL) {
        delete tensor_iterator;
    }
}

double& Tensor::operator[](int* loc) {
    int index = 0;
    for (int i = 0; i < order; i++) {
        index += strides[i] * loc[i];
    }
    return values[index];
}

TensorIterator* Tensor::iterator() {
    tensor_iterator->reset();
    return tensor_iterator;
}

}  // namespace fluoroseq
