/******************************************************************************\
* Author: Matthew Beauregard Smith                                             *
* Affiliation: The University of Texas at Austin                               *
* Department: Oden Institute and Institute for Cellular and Molecular Biology  *
* PI: Edward Marcotte                                                          *
* Project: Protein Fluorosequencing                                            *
\******************************************************************************/

// Defining symbols from header:
#include "binomial-transition.h"

// Local project headers:
#include "tensor/tensor.h"
#include "tensor/vector.h"

namespace fluoroseq {

BinomialTransition::BinomialTransition(double q) : q(q) {
    length = 1;
    size = 1;
    values.resize(size);
    values[0] = 1.0;
}

void BinomialTransition::reserve(int max_n) {
    if (max_n + 1 <= length) {
        return;
    }
    int prev_length = length;
    length = max_n + 1;
    size = length * (length + 1) / 2;
    values.resize(size);
    double p = (double)1 - q;
    for (int i = prev_length; i < length; i++) {
        prob(i, 0) = prob(i - 1, 0) * q;
        for (int j = 1; j < i; j++) {
            prob(i, j) = prob(i - 1, j) * q + prob(i - 1, j - 1) * p;
        }
        prob(i, i) = prob(i - 1, i - 1) * p;
    }
}

double& BinomialTransition::prob(int from, int to) {
    return values[from * (from + 1) / 2 + to];
}

double BinomialTransition::prob(int from, int to) const {
    return values[from * (from + 1) / 2 + to];
}

void BinomialTransition::forward(const Tensor& input,
                                 int channel,
                                 int edmans,
                                 Tensor* output) const {
    int vector_stride = input.strides[1 + channel];
    int vector_length = input.shape[1 + channel];
    int outer_stride = vector_stride * vector_length;
    int outer_max = input.strides[0] * (edmans + 1);
    for (int outer = 0; outer < outer_max; outer += outer_stride) {
        for (int inner = 0; inner < vector_stride; inner++) {
            const Vector inv(
                    vector_length, vector_stride, &input.values[outer + inner]);
            Vector outv(vector_length,
                        vector_stride,
                        &output->values[outer + inner]);
            this->forward(inv, &outv);
        }
    }
}

void BinomialTransition::forward(const Vector& input, Vector* output) const {
    for (int to = 0; to < output->length; to++) {
        double v_to = 0.0;
        for (int from = to; from < input.length; from++) {
            v_to += prob(from, to) * input[from];
        }
        (*output)[to] = v_to;
    }
}

}  // namespace fluoroseq
