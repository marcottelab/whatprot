// Author: Matthew Beauregard Smith (UT Austin)
#include "vector.h"

namespace fluoroseq {

Vector::Vector(int length, int stride, double* values)
        : length(length), stride(stride), values(values) {}

double& Vector::operator[](int i) {
    return values[i * stride];
}

}  // namespace fluoroseq
