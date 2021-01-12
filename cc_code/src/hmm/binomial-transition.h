/******************************************************************************\
* Author: Matthew Beauregard Smith                                             *
* Affiliation: The University of Texas at Austin                               *
* Department: Oden Institute and Institute for Cellular and Molecular Biology  *
* PI: Edward Marcotte                                                          *
* Project: Protein Fluorosequencing                                            *
\******************************************************************************/

#ifndef FLUOROSEQ_HMM_BINOMIAL_TRANSITION
#define FLUOROSEQ_HMM_BINOMIAL_TRANSITION

// Standard C++ library headers:
#include <vector>

// Local project headers:
#include "tensor/tensor.h"
#include "tensor/vector.h"

namespace fluoroseq {

class BinomialTransition {
public:
    BinomialTransition(double q);
    void reserve(int max_n);
    double& prob(int from, int to);
    double prob(int from, int to) const;
    void forward(const Tensor& input, int channel, int edmans, Tensor* output) const;
    void forward(const Vector& input, Vector* output) const;

    std::vector<double> values;
    const double q;
    int length;  // length of array in one dimension.
    int size;  // length of values.
};

}  // namespace fluoroseq

#endif  // FLUOROSEQ_HMM_BINOMIAL_TRANSITION