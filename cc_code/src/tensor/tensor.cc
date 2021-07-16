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
#include <initializer_list>

// Local project headers:
#include "tensor/const-tensor-iterator.h"
#include "tensor/const-tensor-vector-iterator.h"
#include "tensor/tensor-iterator.h"
#include "tensor/tensor-vector-iterator.h"
#include "util/kd-range.h"

namespace whatprot {

namespace {
using std::copy;
using std::initializer_list;
}  // namespace

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

double& Tensor::operator[](const unsigned int* loc) {
    unsigned int index = 0;
    for (unsigned int i = 0; i < order; i++) {
        index += strides[i] * loc[i];
    }
    return values[index];
}

double& Tensor::operator[](initializer_list<unsigned int> loc) {
    return (*this)[loc.begin()];
}

TensorIterator* Tensor::iterator(const KDRange& range) {
    return new TensorIterator(order, range, shape, size, values);
}

ConstTensorIterator* Tensor::const_iterator(const KDRange& range) const {
    return new ConstTensorIterator(order, range, shape, size, values);
}

TensorVectorIterator* Tensor::vector_iterator(const KDRange& range,
                                              unsigned int vector_dimension) {
    return new TensorVectorIterator(
            order, range, shape, strides, size, values, vector_dimension);
}

ConstTensorVectorIterator* Tensor::const_vector_iterator(
        const KDRange& range, unsigned int vector_dimension) const {
    return new ConstTensorVectorIterator(
            order, range, shape, strides, size, values, vector_dimension);
}

double Tensor::sum() const {
    double total = 0.0;
    for (unsigned int i = 0; i < size; i++) {
        total += values[i];
    }
    return total;
}

}  // namespace whatprot
