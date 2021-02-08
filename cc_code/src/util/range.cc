/******************************************************************************\
* Author: Matthew Beauregard Smith                                             *
* Affiliation: The University of Texas at Austin                               *
* Department: Oden Institute and Institute for Cellular and Molecular Biology  *
* PI: Edward Marcotte                                                          *
* Project: Protein Fluorosequencing                                            *
\******************************************************************************/

// Defining symbols from header:
#include "range.h"

namespace whatprot {

Range::Range(int max) : min(0), max(max) {}

Range::Range(int min, int max) : min(min), max(max) {}

RangeIterator Range::begin() {
    return RangeIterator(min);
}

RangeIterator Range::end() {
    return RangeIterator(max);
}

RangeIterator::RangeIterator(int index) : index(index) {}

void RangeIterator::operator++() {
    index++;
}

bool RangeIterator::operator!=(const RangeIterator& other) {
    return (index != other.index);
}

int RangeIterator::operator*() {
    return index;
}

}  // namespace whatprot
