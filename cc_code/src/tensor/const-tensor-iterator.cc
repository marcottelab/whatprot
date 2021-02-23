/******************************************************************************\
* Author: Matthew Beauregard Smith                                             *
* Affiliation: The University of Texas at Austin                               *
* Department: Oden Institute and Institute for Cellular and Molecular Biology  *
* PI: Edward Marcotte                                                          *
* Project: Protein Fluorosequencing                                            *
\******************************************************************************/

// Defining symbols from header:
#include "const-tensor-iterator.h"

namespace whatprot {

ConstTensorIterator::ConstTensorIterator(int order,
                                         const int* shape,
                                         int size,
                                         const double* values)
        : order(order), shape(shape), size(size), values(values), index(0) {
    loc = new int[order]();
    reset();
}

ConstTensorIterator::~ConstTensorIterator() {
    delete[] loc;
}

void ConstTensorIterator::reset() {
    for (int o = 0; o < order; o++) {
        loc[o] = 0;
    }
    index = 0;
}

void ConstTensorIterator::advance() {
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

double ConstTensorIterator::get() const {
    return values[index];
}

bool ConstTensorIterator::done() {
    return index == size;
}

}  // namespace whatprot
