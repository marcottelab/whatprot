/******************************************************************************\
* Author: Matthew Beauregard Smith                                             *
* Affiliation: The University of Texas at Austin                               *
* Department: Oden Institute and Institute for Cellular and Molecular Biology  *
* PI: Edward Marcotte                                                          *
* Project: Protein Fluorosequencing                                            *
\******************************************************************************/

#ifndef FLUOROSEQ_HMM_BINOMIAL_TRANSITION_H
#define FLUOROSEQ_HMM_BINOMIAL_TRANSITION_H

// Standard C++ library headers:
#include <vector>

// Local project headers:
#include "hmm/step.h"
#include "tensor/tensor.h"
#include "tensor/vector.h"

namespace fluoroseq {

class BinomialTransition : public Step {
public:
    BinomialTransition(double q, int channel);
    void reserve(int max_n);
    double& prob(int from, int to);
    double prob(int from, int to) const;
    virtual void forward(const Tensor& input,
                         int* edmans,
                         Tensor* output) const override;
    void forward(const Vector& input, Vector* output) const;
    virtual void backward(const Tensor& input,
                          int* edmans,
                          Tensor* output) const override;
    void backward(const Vector& input, Vector* output) const;

    std::vector<double> values;
    const double q;
    int channel;
    int length;  // length of array in one dimension.
    int size;  // length of values.
};

}  // namespace fluoroseq

#endif  // FLUOROSEQ_HMM_BINOMIAL_TRANSITION_H