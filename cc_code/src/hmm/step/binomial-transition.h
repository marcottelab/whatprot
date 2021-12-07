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
#include "hmm/state-vector/peptide-state-vector.h"
#include "hmm/step/peptide-step.h"
#include "parameterization/fit/parameter-fitter.h"
#include "tensor/tensor.h"
#include "tensor/vector.h"
#include "util/kd-range.h"

namespace whatprot {

class BinomialTransition : public PeptideStep {
public:
    BinomialTransition(double q, int channel);
    void reserve(unsigned int max_n);
    double& prob(unsigned int from, unsigned int to);
    double prob(unsigned int from, unsigned int to) const;
    virtual void prune_forward(KDRange* range, bool* allow_detached) override;
    virtual void prune_backward(KDRange* range, bool* allow_detached) override;
    virtual void forward(unsigned int* num_edmans,
                         PeptideStateVector* psv) const override;
    void forward(Vector* v) const;
    virtual void backward(const PeptideStateVector& input,
                          unsigned int* num_edmans,
                          PeptideStateVector* output) const override;
    void backward(const Vector& input, Vector* output) const;
    void improve_fit(const PeptideStateVector& forward_psv,
                     const PeptideStateVector& backward_psv,
                     const PeptideStateVector& next_backward_psv,
                     unsigned int num_edmans,
                     double probability,
                     ParameterFitter* fitter) const;
    void improve_fit(const Vector& forward_vector,
                     const Vector& backward_vector,
                     const Vector& next_backward_vector,
                     double probability,
                     ParameterFitter* fitter) const;
    KDRange forward_range;
    KDRange backward_range;
    std::vector<double> values;
    const double q;
    int channel;
    unsigned int length;  // length of array in one dimension.
    unsigned int size;  // length of values.
};

}  // namespace whatprot

#endif  // WHATPROT_HMM_STEP_BINOMIAL_TRANSITION_H