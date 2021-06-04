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
#include "hmm/state-vector/peptide-state-vector.h"
#include "parameterization/fit/parameter-fitter.h"
#include "tensor/tensor.h"
#include "tensor/vector.h"

namespace whatprot {

BinomialTransition::BinomialTransition(double q, int channel)
        : q(q), channel(channel) {
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

void BinomialTransition::forward(int* num_edmans,
                                 PeptideStateVector* psv) const {
    int vector_stride = psv->tensor.strides[1 + channel];
    int vector_length = psv->tensor.shape[1 + channel];
    int outer_stride = vector_stride * vector_length;
    int outer_max = psv->tensor.strides[0] * (*num_edmans + 1);
    for (int outer = 0; outer < outer_max; outer += outer_stride) {
        for (int inner = 0; inner < vector_stride; inner++) {
            Vector v(vector_length,
                     vector_stride,
                     &psv->tensor.values[outer + inner]);
            this->forward(&v);
        }
    }
}

void BinomialTransition::forward(Vector* v) const {
    for (int to = 0; to < v->length; to++) {
        double v_to = 0.0;
        for (int from = to; from < v->length; from++) {
            v_to += prob(from, to) * (*v)[from];
        }
        (*v)[to] = v_to;
    }
}

void BinomialTransition::backward(const PeptideStateVector& input,
                                  int* num_edmans,
                                  PeptideStateVector* output) const {
    int vector_stride = input.tensor.strides[1 + channel];
    int vector_length = input.tensor.shape[1 + channel];
    int outer_stride = vector_stride * vector_length;
    int outer_max = input.tensor.strides[0] * (*num_edmans + 1);
    for (int outer = 0; outer < outer_max; outer += outer_stride) {
        for (int inner = 0; inner < vector_stride; inner++) {
            const Vector inv(vector_length,
                             vector_stride,
                             &input.tensor.values[outer + inner]);
            Vector outv(vector_length,
                        vector_stride,
                        &output->tensor.values[outer + inner]);
            this->backward(inv, &outv);
        }
    }
}

void BinomialTransition::backward(const Vector& input, Vector* output) const {
    for (int from = output->length - 1; from >= 0; from--) {
        double v_from = 0.0;
        for (int to = 0; to <= from; to++) {
            v_from += prob(from, to) * input[to];
        }
        (*output)[from] = v_from;
    }
}

void BinomialTransition::improve_fit(
        const PeptideStateVector& forward_psv,
        const PeptideStateVector& backward_psv,
        const PeptideStateVector& next_backward_psv,
        int num_edmans,
        double probability,
        ParameterFitter* fitter) const {
    int vector_stride = forward_psv.tensor.strides[1 + channel];
    int vector_length = forward_psv.tensor.shape[1 + channel];
    int outer_stride = vector_stride * vector_length;
    int outer_max = forward_psv.tensor.strides[0] * (num_edmans + 1);
    for (int outer = 0; outer < outer_max; outer += outer_stride) {
        for (int inner = 0; inner < vector_stride; inner++) {
            const Vector fv(vector_length,
                            vector_stride,
                            &forward_psv.tensor.values[outer + inner]);
            const Vector bv(vector_length,
                            vector_stride,
                            &backward_psv.tensor.values[outer + inner]);
            const Vector nbv(vector_length,
                             vector_stride,
                             &next_backward_psv.tensor.values[outer + inner]);
            this->improve_fit(fv, bv, nbv, probability, fitter);
        }
    }
}

void BinomialTransition::improve_fit(const Vector& forward_vector,
                                     const Vector& backward_vector,
                                     const Vector& next_backward_vector,
                                     double probability,
                                     ParameterFitter* fitter) const {
    // Note that we can ignore when starting location (from) is 0 because then
    // there are no dyes, it's irrelevant. We would be adding 0s.
    for (int from = forward_vector.length - 1; from > 0; from--) {
        double p_state =
                forward_vector[from] * backward_vector[from] / probability;
        fitter->denominator += p_state * (double)from;
        // We can ignore when to and from are equal, because no dyes are lost
        // then, so it gives us nothing else for the numerator; we would be
        // adding zero.
        for (int to = 0; to < from; to++) {
            double p_transition = forward_vector[from] * prob(from, to)
                                  * next_backward_vector[to] / probability;
            fitter->numerator += p_transition * (double)(from - to);
        }
    }
}

}  // namespace whatprot
