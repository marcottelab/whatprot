/******************************************************************************\
* Author: Matthew Beauregard Smith                                             *
* Affiliation: The University of Texas at Austin                               *
* Department: Oden Institute and Institute for Cellular and Molecular Biology  *
* PI: Edward Marcotte                                                          *
* Project: Protein Fluorosequencing                                            *
\******************************************************************************/

#ifndef FLUOROSEQ_FWD_ALG_BINOMIAL_TRANSITION
#define FLUOROSEQ_FWD_ALG_BINOMIAL_TRANSITION

// Standard C++ library headers:
#include <vector>

// Local project headers:
#include "tensor/tensor.h"
#include "tensor/vector.h"

namespace fluoroseq {

class BinomialTransition {
  public:
    BinomialTransition(double q);
    void reserve(int max_n) const;
    double& prob(int from, int to) const;
    void operator()(Tensor* tensor, int channel, int edmans) const;
    void operator()(Vector* v) const;

    mutable std::vector<double> values;
    const double q;
    mutable int length;  // length of array in one dimension.
    mutable int size;  // length of values.
};

}  // namespace fluoroseq

#endif  // FLUOROSEQ_FWD_ALG_BINOMIAL_TRANSITION