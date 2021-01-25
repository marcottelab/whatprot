/******************************************************************************\
* Author: Matthew Beauregard Smith                                             *
* Affiliation: The University of Texas at Austin                               *
* Department: Oden Institute and Institute for Cellular and Molecular Biology  *
* PI: Edward Marcotte                                                          *
* Project: Protein Fluorosequencing                                            *
\******************************************************************************/

#ifndef FLUOROSEQ_HMM_START_H
#define FLUOROSEQ_HMM_START_H

// Local project headers:
#include "hmm/step.h"
#include "tensor/tensor.h"

namespace fluoroseq {

class Start : public Step {
public:
    virtual void forward(const Tensor& input,
                         int* edmans,
                         Tensor* output) const override;
    virtual void backward(const Tensor& input,
                          int* edmans,
                          Tensor* output) const override;
};

}  // namespace fluoroseq

#endif  // FLUOROSEQ_HMM_START_H