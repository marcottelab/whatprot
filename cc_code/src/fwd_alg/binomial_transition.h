// Author: Matthew Beauregard Smith
#ifndef FLUOROSEQ_FWD_ALG_BINOMIAL_TRANSITION
#define FLUOROSEQ_FWD_ALG_BINOMIAL_TRANSITION

#include "tensor/tensor.h"
#include "tensor/vector.h"

namespace fluoroseq {

class BinomialTransition {
public:
    BinomialTransition(int max_n, double q, int max_failed_edmans);
    ~BinomialTransition();
    double& prob(int from, int to);
    double prob(int from, int to) const;
    void operator()(Tensor* tensor, int channel, int edmans) const;
    void operator()(Vector* v) const;

    double* values;
    int length;  // length of array in one dimension.
    int size;  // length of values.
    int max_failed_edmans;
};

}  // namespace fluoroseq

#endif  // FLUOROSEQ_FWD_ALG_BINOMIAL_TRANSITION