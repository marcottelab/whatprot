/******************************************************************************\
* Author: Matthew Beauregard Smith                                             *
* Affiliation: The University of Texas at Austin                               *
* Department: Oden Institute and Institute for Cellular and Molecular Biology  *
* PI: Edward Marcotte                                                          *
* Project: Protein Fluorosequencing                                            *
\******************************************************************************/

// Defining symbols from header:
#include "binomial-transition.h"

// Standard C++ library headers:
#include <limits>

// Local project headers:
#include "hmm/state-vector/peptide-state-vector.h"
#include "parameterization/fit/parameter-fitter.h"
#include "tensor/tensor.h"
#include "tensor/vector.h"
#include "util/kd-range.h"

namespace whatprot {

namespace {
using std::numeric_limits;
}  // namespace

BinomialTransition::BinomialTransition(double q, int channel)
        : q(q), channel(channel) {
    length = 1;
    size = 1;
    values.resize(size);
    values[0] = 1.0;
}

void BinomialTransition::reserve(unsigned int max_n) {
    if (max_n + 1 <= length) {
        return;
    }
    unsigned int prev_length = length;
    length = max_n + 1;
    size = length * (length + 1) / 2;
    values.resize(size);
    double p = (double)1 - q;
    for (unsigned int i = prev_length; i < length; i++) {
        prob(i, 0) = prob(i - 1, 0) * q;
        for (unsigned int j = 1; j < i; j++) {
            prob(i, j) = prob(i - 1, j) * q + prob(i - 1, j - 1) * p;
        }
        prob(i, i) = prob(i - 1, i - 1) * p;
    }
}

double& BinomialTransition::prob(unsigned int from, unsigned int to) {
    return values[from * (from + 1) / 2 + to];
}

double BinomialTransition::prob(unsigned int from, unsigned int to) const {
    return values[from * (from + 1) / 2 + to];
}

void BinomialTransition::prune_forward(KDRange* range, bool* allow_detached) {
    forward_range = *range;
    range->min[1 + channel] = 0;
    backward_range = *range;
}

void BinomialTransition::prune_backward(KDRange* range, bool* allow_detached) {
    backward_range = backward_range.intersect(*range);
    *range = backward_range;
    range->max[1 + channel] = numeric_limits<unsigned int>::max();
    forward_range = forward_range.intersect(*range);
    *range = forward_range;
}

void BinomialTransition::forward(const PeptideStateVector& input,
                                 unsigned int* num_edmans,
                                 PeptideStateVector* output) const {
    // Mismatched range is OK, range should only differ on the channel being
    // processed.
    ConstTensorVectorIterator* in_itr =
            input.tensor.const_vector_iterator(forward_range, 1 + channel);
    TensorVectorIterator* out_itr =
            output->tensor.vector_iterator(backward_range, 1 + channel);
    while (!in_itr->done()) {
        const Vector* in_v = in_itr->get();
        Vector* out_v = out_itr->get();
        this->forward(*in_v, out_v);
        delete in_v;
        delete out_v;
        in_itr->advance();
        out_itr->advance();
    }
    delete in_itr;
    delete out_itr;
    output->range = backward_range;
}

void BinomialTransition::forward(const Vector& input, Vector* output) const {
    unsigned int to_min = backward_range.min[1 + channel];
    unsigned int to_max = backward_range.max[1 + channel];
    for (unsigned int to = to_min; to < to_max; to++) {
        double v_to = 0.0;
        unsigned int from_min = std::max(to, forward_range.min[1 + channel]);
        unsigned int from_max = forward_range.max[1 + channel];
        for (unsigned int from = from_min; from < from_max; from++) {
            v_to += prob(from, to) * input[from];
        }
        (*output)[to] = v_to;
    }
}

void BinomialTransition::backward(const PeptideStateVector& input,
                                  unsigned int* num_edmans,
                                  PeptideStateVector* output) const {
    ConstTensorVectorIterator* in_itr =
            input.tensor.const_vector_iterator(backward_range, 1 + channel);
    TensorVectorIterator* out_itr =
            output->tensor.vector_iterator(forward_range, 1 + channel);
    while (!out_itr->done()) {
        const Vector* in_v = in_itr->get();
        Vector* out_v = out_itr->get();
        this->backward(*in_v, out_v);
        delete in_v;
        delete out_v;
        in_itr->advance();
        out_itr->advance();
    }
    delete in_itr;
    delete out_itr;
    output->range = forward_range;
}

void BinomialTransition::backward(const Vector& input, Vector* output) const {
    int from_min = forward_range.min[1 + channel];
    int from_max = forward_range.max[1 + channel];
    for (int from = from_max - 1; from >= from_min; from--) {
        double v_from = 0.0;
        int to_min = backward_range.min[1 + channel];
        int to_max = std::min(from, (int)backward_range.max[1 + channel] - 1);
        for (int to = to_min; to <= to_max; to++) {
            v_from += prob(from, to) * input[to];
        }
        (*output)[from] = v_from;
    }
}

void BinomialTransition::improve_fit(
        const PeptideStateVector& forward_psv,
        const PeptideStateVector& backward_psv,
        const PeptideStateVector& next_backward_psv,
        unsigned int num_edmans,
        double probability,
        ParameterFitter* fitter) const {
    ConstTensorVectorIterator* f_itr = forward_psv.tensor.const_vector_iterator(
            forward_range, 1 + channel);
    ConstTensorVectorIterator* b_itr =
            backward_psv.tensor.const_vector_iterator(forward_range,
                                                      1 + channel);
    ConstTensorVectorIterator* nb_itr =
            next_backward_psv.tensor.const_vector_iterator(backward_range,
                                                           1 + channel);
    while (!f_itr->done()) {
        const Vector* f_v = f_itr->get();
        const Vector* b_v = b_itr->get();
        const Vector* nb_v = nb_itr->get();
        this->improve_fit(*f_v, *b_v, *nb_v, probability, fitter);
        delete f_v;
        delete b_v;
        delete nb_v;
        f_itr->advance();
        b_itr->advance();
        nb_itr->advance();
    }
    delete f_itr;
    delete b_itr;
    delete nb_itr;
}

void BinomialTransition::improve_fit(const Vector& forward_vector,
                                     const Vector& backward_vector,
                                     const Vector& next_backward_vector,
                                     double probability,
                                     ParameterFitter* fitter) const {
    int from_min = forward_range.min[1 + channel];
    int from_max = forward_range.max[1 + channel];
    // Note that we can ignore when starting location (from) is 0 because then
    // there are no dyes, it's irrelevant. We would be adding 0s.
    for (int from = from_max - 1; from >= std::max(from_min, 1); from--) {
        double p_state =
                forward_vector[from] * backward_vector[from] / probability;
        fitter->denominator += p_state * (double)from;
        int to_min = backward_range.min[1 + channel];
        int to_max = backward_range.max[1 + channel];
        // We can ignore when to and from are equal, because no dyes are lost
        // then, so it gives us nothing else for the numerator; we would be
        // adding zero.
        for (int to = to_min; to < std::min(to_max, from); to++) {
            double p_transition = forward_vector[from] * prob(from, to)
                                  * next_backward_vector[to] / probability;
            fitter->numerator += p_transition * (double)(from - to);
        }
    }
}

}  // namespace whatprot
