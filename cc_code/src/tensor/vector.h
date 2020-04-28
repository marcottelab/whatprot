// Author: Matthew Beauregard Smith (UT Austin)
#ifndef FLUOROSEQ_TENSOR_VECTOR_H
#define FLUOROSEQ_TENSOR_VECTOR_H

namespace fluoroseq {

class Vector {
public:
    Vector(int length, int stride, double* values);
    double& operator[](int i);

    double* values;
    int length;
    int stride;
};

}  // namespace fluoroseq

#endif  // FLUOROSEQ_TENSOR_VECTOR_H
