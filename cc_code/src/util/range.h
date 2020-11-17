/******************************************************************************\
* Author: Matthew Beauregard Smith                                             *
* Affiliation: The University of Texas at Austin                               *
* Department: Oden Institute and Institute for Cellular and Molecular Biology  *
* PI: Edward Marcotte                                                          *
* Project: Protein Fluorosequencing                                            *
\******************************************************************************/

#ifndef FLUOROSEQ_UTIL_RANGE_H
#define FLUOROSEQ_UTIL_RANGE_H

namespace fluoroseq {

class RangeIterator;

class Range {
  public:
    Range(int max);
    Range(int min, int max);
    RangeIterator begin();
    RangeIterator end();

    int min;
    int max;
};

class RangeIterator {
  public:
    RangeIterator(int index);
    void operator++();
    bool operator!=(const RangeIterator& other);
    int operator*();

    int index;
};

}  // namespace fluoroseq

#endif  // FLUOROSEQ_UTIL_RANGE_H
