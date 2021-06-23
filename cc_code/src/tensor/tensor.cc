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
#include "tensor/const-tensor-iterator.h"
#include "tensor/tensor-iterator.h"

namespace whatprot {

namespace {
using std::copy;
}

Tensor::Tensor(unsigned int order, const unsigned int* shape) : order(order) {
    this->shape = new unsigned int[order];
    copy(shape, shape + order, this->shape);
    size = 1;
    strides = new int[order];
    for (int i = order - 1; i >= 0; i--) {
        strides[i] = size;
        size *= shape[i];
    }
    values = new double[size]();
}

Tensor::Tensor(Tensor&& other)
        : values(other.values),
          shape(other.shape),
          strides(other.strides),
          size(other.size),
          order(other.order) {
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
}

double& Tensor::operator[](unsigned int* loc) {
    unsigned int index = 0;
    for (unsigned int i = 0; i < order; i++) {
        index += strides[i] * loc[i];
    }
    return values[index];
}

TensorIterator* Tensor::iterator(const unsigned int* min,
                                 const unsigned int* max) {
    return new TensorIterator(order, min, max, shape, size, values);
}

ConstTensorIterator* Tensor::const_iterator(const unsigned int* min,
                                            const unsigned int* max) const {
    return new ConstTensorIterator(order, min, max, shape, size, values);
}

double Tensor::sum() const {
    double total = 0.0;
    for (unsigned int i = 0; i < size; i++) {
        total += values[i];
    }
    return total;
}

}  // namespace whatprot
