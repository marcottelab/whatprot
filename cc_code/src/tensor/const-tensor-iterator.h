/******************************************************************************\
* Author: Matthew Beauregard Smith                                             *
* Affiliation: The University of Texas at Austin                               *
* Department: Oden Institute and Institute for Cellular and Molecular Biology  *
* PI: Edward Marcotte                                                          *
* Project: Protein Fluorosequencing                                            *
\******************************************************************************/

#ifndef FLUOROSEQ_TENSOR_CONST_TENSOR_ITERATOR_H
#define FLUOROSEQ_TENSOR_CONST_TENSOR_ITERATOR_H

namespace fluoroseq {

class ConstTensorIterator {
public:
    ConstTensorIterator(int order, const int* shape, int size, const double* values);
    ~ConstTensorIterator();
    void reset();
    void advance();
    double get() const;
    bool done();

    const double* values;  // not owned
    const int* shape;  // not owned
    int* loc;
    const int order;
    int index;  // current index directly into values
    const int size;  // length of values
};

}  // namespace fluoroseq

#endif  // FLUOROSEQ_TENSOR_CONST_TENSOR_ITERATOR_H