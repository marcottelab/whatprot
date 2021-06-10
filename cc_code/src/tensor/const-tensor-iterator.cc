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

ConstTensorIterator::ConstTensorIterator(unsigned int order,
                                         const unsigned int* shape,
                                         unsigned int size,
                                         const double* values)
        : values(values), shape(shape), order(order), index(0), size(size) {
    loc = new unsigned int[order]();
    reset();
}

ConstTensorIterator::~ConstTensorIterator() {
    delete[] loc;
}

void ConstTensorIterator::reset() {
    for (unsigned int o = 0; o < order; o++) {
        loc[o] = 0;
    }
    index = 0;
}

void ConstTensorIterator::advance() {
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

double ConstTensorIterator::get() const {
    return values[index];
}

bool ConstTensorIterator::done() {
    return index == size;
}

}  // namespace whatprot
