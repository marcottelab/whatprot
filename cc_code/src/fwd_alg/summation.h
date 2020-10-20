// Author: Matthew Beauregard Smith (UT Austin)
#ifndef FLUOROSEQ_FWD_ALG_SUMMATION_H
#define FLUOROSEQ_FWD_ALG_SUMMATION_H

#include "tensor/tensor.h"

namespace fluoroseq {

class Summation {
public:
    Summation();
    double operator()(Tensor* tensor, int timestep) const;
};

}  // namespace fluoroseq

#endif  // FLUOROSEQ_FWD_ALG_SUMMATION_H