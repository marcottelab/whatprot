// Author: Matthew Beauregard Smith (UT Austin)
#include "tensor_iterator.h"

namespace fluoroseq {

TensorIterator::TensorIterator(int order, int* shape, int size, double* values)
        : order(order), shape(shape), size(size), values(values), index(0) {
    loc = new int[order]();
}

TensorIterator::~TensorIterator() {
    delete[] loc;
}

void TensorIterator::reset() {
    for (int o = 0; o < order; o++) {
        loc[o] = 0;
    }
    index = 0;
}

void TensorIterator::advance() {
    index++;
    for (int o = order - 1; o >= 0; o--) {
        loc[o]++;
        if (loc[o] < shape[o]) {
            break;
        } else {
            loc[o] = 0;
        }
    }
}

double* TensorIterator::get() {
    return &values[index];
}

bool TensorIterator::done() {
    return index == size;
}

}  // namespace fluoroseq
