/******************************************************************************\
* Author: Matthew Beauregard Smith                                             *
* Affiliation: The University of Texas at Austin                               *
* Department: Oden Institute and Institute for Cellular and Molecular Biology  *
* PI: Edward Marcotte                                                          *
* Project: Protein Fluorosequencing                                            *
\******************************************************************************/

// Defining symbols from header:
#include "tensor-iterator.h"

namespace whatprot {

TensorIterator::TensorIterator(unsigned int order,
                               unsigned int* shape,
                               unsigned int size,
                               double* values)
        : values(values), shape(shape), order(order), index(0), size(size) {
    loc = new unsigned int[order]();
    reset();
}

TensorIterator::~TensorIterator() {
    delete[] loc;
}

void TensorIterator::reset() {
    for (unsigned int o = 0; o < order; o++) {
        loc[o] = 0;
    }
    index = 0;
}

void TensorIterator::advance() {
    index++;
    // Need to use signed int o to detect when out of entries to decrease.
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

}  // namespace whatprot
