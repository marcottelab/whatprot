/******************************************************************************\
* Author: Matthew Beauregard Smith                                             *
* Affiliation: The University of Texas at Austin                               *
* Department: Oden Institute and Institute for Cellular and Molecular Biology  *
* PI: Edward Marcotte                                                          *
* Project: Protein Fluorosequencing                                            *
\******************************************************************************/

#ifndef WHATPROT_TENSOR_CONST_TENSOR_ITERATOR_H
#define WHATPROT_TENSOR_CONST_TENSOR_ITERATOR_H

namespace whatprot {

class ConstTensorIterator {
public:
    ConstTensorIterator(unsigned int order,
                        const unsigned int* shape,
                        unsigned int size,
                        const double* values);
    ~ConstTensorIterator();
    void reset();
    void advance();
    double get() const;
    bool done();

    const double* values;  // not owned
    const unsigned int* shape;  // not owned
    unsigned int* loc;
    const unsigned int order;
    unsigned int index;  // current index directly into values
    const unsigned int size;  // length of values
};

}  // namespace whatprot

#endif  // WHATPROT_TENSOR_CONST_TENSOR_ITERATOR_H