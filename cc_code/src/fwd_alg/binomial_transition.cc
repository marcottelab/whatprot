// Author: Matthew Beauregard Smith
#include "binomial_transition.h"

#include <algorithm>

#include "tensor/tensor.h"
#include "tensor/vector.h"

namespace fluoroseq {

namespace {
using std::max;
}  // namespace

BinomialTransition::BinomialTransition(double q, int max_failed_edmans)
            : q(q), max_failed_edmans(max_failed_edmans) {
    length = 1;
    size = 1;
    values.reserve(size);
    values[0] = 1.0;
}

void BinomialTransition::reserve(int max_n) const {
    int prev_length = length;
    length = max_n + 1;
    size = length * (length + 1) / 2;
    values.reserve(size);
    double p = (double) 1 - q;
    for (int i = prev_length; i < length; i++) {
        prob(i, 0) = prob(i - 1, 0) * q;
        for (int j = 1; j < i; j++) {
            prob(i, j) = prob(i - 1, j) * q + prob(i - 1, j - 1) * p;
        }
        prob(i, i) = prob(i - 1, i - 1) * p;
    }
}

double& BinomialTransition::prob(int from, int to) const {
    return values[from * (from + 1) / 2 + to];
}

void BinomialTransition::operator()(Tensor* tensor,
                                    int channel,
                                    int edmans) const {
    int vector_stride = tensor->strides[1 + channel];
    int vector_length = tensor->shape[1 + channel];
    reserve(vector_length);
    int outer_stride = vector_stride * vector_length;
    int outer_min = max(0, tensor->strides[0] * (edmans - max_failed_edmans));
    int outer_max = tensor->strides[0] * (edmans + 1);
    for (int outer = outer_min; outer < outer_max; outer += outer_stride) {
        for (int inner = 0; inner < vector_stride; inner++) {
            Vector v(vector_length,
                     vector_stride,
                     &tensor->values[outer + inner]);
            (*this)(&v);
        }
    }
}

void BinomialTransition::operator()(Vector* v) const {
    for (int to = 0; to < v->length; to++) {
        double v_to = 0.0;
        for (int from = to; from < v->length; from++) {
            v_to += prob(from, to) * (*v)[from];
        }
        (*v)[to] = v_to;
    }
}

}  // namespace fluoroseq
