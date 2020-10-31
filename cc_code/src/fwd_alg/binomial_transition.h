// Author: Matthew Beauregard Smith
#ifndef FLUOROSEQ_FWD_ALG_BINOMIAL_TRANSITION
#define FLUOROSEQ_FWD_ALG_BINOMIAL_TRANSITION

#include <vector>

#include "tensor/tensor.h"
#include "tensor/vector.h"

namespace fluoroseq {

class BinomialTransition {
public:
    BinomialTransition(double q);
    void reserve(int max_n) const;
    void reserve(int max_n);
    double prob(int from, int to) const;
    double& prob(int from, int to);
    void operator()(Tensor* tensor, int channel, int edmans) const;
    void operator()(Vector* v) const;

    mutable std::vector<double> values;
    const double q;
    mutable int length;  // length of array in one dimension.
    mutable int size;  // length of values.
};

}  // namespace fluoroseq

#endif  // FLUOROSEQ_FWD_ALG_BINOMIAL_TRANSITION