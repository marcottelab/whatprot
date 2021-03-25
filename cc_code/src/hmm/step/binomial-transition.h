/******************************************************************************\
* Author: Matthew Beauregard Smith                                             *
* Affiliation: The University of Texas at Austin                               *
* Department: Oden Institute and Institute for Cellular and Molecular Biology  *
* PI: Edward Marcotte                                                          *
* Project: Protein Fluorosequencing                                            *
\******************************************************************************/

#ifndef WHATPROT_HMM_STEP_BINOMIAL_TRANSITION_H
#define WHATPROT_HMM_STEP_BINOMIAL_TRANSITION_H

// Standard C++ library headers:
#include <vector>

// Local project headers:
#include "hmm/fit/parameter-fitter.h"
#include "hmm/state-vector/peptide-state-vector.h"
#include "hmm/step/step.h"
#include "tensor/tensor.h"
#include "tensor/vector.h"

namespace whatprot {

class BinomialTransition : public Step<PeptideStateVector> {
public:
    BinomialTransition(double q, int channel);
    void reserve(int max_n);
    double& prob(int from, int to);
    double prob(int from, int to) const;
    virtual void forward(int* num_edmans, PeptideStateVector* psv) const override;
    void forward(Vector* v) const;
    virtual void backward(const PeptideStateVector& input, int* num_edmans,
                          PeptideStateVector* output) const override;
    void backward(const Vector& input, Vector* output) const;
    void improve_fit(const PeptideStateVector& forward_psv,
                     const PeptideStateVector& backward_psv,
                     const PeptideStateVector& next_backward_psv, int num_edmans,
                     double probability,
                     ParameterFitter* fitter) const;
    void improve_fit(const Vector& forward_vector,
                     const Vector& backward_vector,
                     const Vector& next_backward_vector,
                     double probability,
                     ParameterFitter* fitter) const;
    std::vector<double> values;
    const double q;
    int channel;
    int length;  // length of array in one dimension.
    int size;  // length of values.
};

}  // namespace whatprot

#endif  // WHATPROT_HMM_STEP_BINOMIAL_TRANSITION_H