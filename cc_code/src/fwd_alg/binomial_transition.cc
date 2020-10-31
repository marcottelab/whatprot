// Author: Matthew Beauregard Smith
#include "binomial_transition.h"

#include "tensor/tensor.h"
#include "tensor/vector.h"

namespace fluoroseq {

BinomialTransition::BinomialTransition(double q) : q(q) {
    length = 1;
    size = 1;
    values.resize(size);
    values[0] = 1.0;
}

void BinomialTransition::reserve(int max_n) const {
    // This is a dirty hack. Generally const_cast<>() should be avoided in C++
    // code. However, we require compatibility with the Intel C++ compiler
    // (icc/icpc). The Intel C++ compiler seems to have a bug when dealing with
    // mutables. Mutable members are the standard way of dealing with cases
    // where the class should be const but some of the data shouldn't. However,
    // trying to return a mutable reference to a mutable member variable from a
    // const function appears to fail with the Intel C++ compiler. It was
    // unfortunately necessary to return a mutable reference from the const
    // prob() function. It was cleaner to hack around this in the reserve()
    // function than in the prob() function. I've left the affected variables
    // marked as mutable, even though this is technically no longer necessary,
    // to make the code clearer.
    BinomialTransition* self = const_cast<BinomialTransition*>(this);
    self->reserve(max_n);
}

void BinomialTransition::reserve(int max_n) {
    if (max_n + 1 <= length) {
        return;
    }
    int prev_length = length;
    length = max_n + 1;
    size = length * (length + 1) / 2;
    values.resize(size);
    double p = (double) 1 - q;
    for (int i = prev_length; i < length; i++) {
        prob(i, 0) = prob(i - 1, 0) * q;
        for (int j = 1; j < i; j++) {
            prob(i, j) = prob(i - 1, j) * q + prob(i - 1, j - 1) * p;
        }
        prob(i, i) = prob(i - 1, i - 1) * p;
    }
}

double BinomialTransition::prob(int from, int to) const {
    return values[from * (from + 1) / 2 + to];
}

double& BinomialTransition::prob(int from, int to) {
    return values[from * (from + 1) / 2 + to];
}

void BinomialTransition::operator()(Tensor* tensor,
                                    int channel,
                                    int edmans) const {
    int vector_stride = tensor->strides[1 + channel];
    int vector_length = tensor->shape[1 + channel];
    reserve(vector_length);
    int outer_stride = vector_stride * vector_length;
    int outer_max = tensor->strides[0] * (edmans + 1);
    for (int outer = 0; outer < outer_max; outer += outer_stride) {
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
