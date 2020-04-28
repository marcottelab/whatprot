// Author: Matthew Beauregard Smith (UT Austin)
#ifndef FLUOROSEQ_TENSOR_TENSOR_ITERATOR_H
#define FLUOROSEQ_TENSOR_TENSOR_ITERATOR_H

namespace fluoroseq {

class TensorIterator {
public:
    TensorIterator(int order, int* shape, int size, double* values);
    ~TensorIterator();
    void advance();
    double* get();
    bool done();

    double* values;  // not owned
    int* shape;  // not owned
    int* loc;
    int order;
    int index;  // current index directly into values
    int size;  // length of values
};

}  // namespace fluoroseq

#endif  // FLUOROSEQ_TENSOR_TENSOR_ITERATOR_H