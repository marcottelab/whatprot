/******************************************************************************\
* Author: Matthew Beauregard Smith                                             *
* Affiliation: The University of Texas at Austin                               *
* Department: Oden Institute and Institute for Cellular and Molecular Biology  *
* PI: Edward Marcotte                                                          *
* Project: Protein Fluorosequencing                                            *
\******************************************************************************/

#ifndef WHATPROT_HMM_STEP_PEPTIDE_STEP_H
#define WHATPROT_HMM_STEP_PEPTIDE_STEP_H

// Local project headers:
#include "hmm/state-vector/peptide-state-vector.h"
#include "hmm/step/step.h"
#include "util/kd-range.h"

namespace whatprot {

// Introduces methods to prune the computation for greater efficiency. This is
// intended to make good use of the near-sparsity of this problem. Caller will
// use prune_forward() before using prune_backward(). Caller will check whether
// the KDRange is empty, the implementations do not need to worry about this.
class PeptideStep : public Step<PeptideStateVector> {
public:
    virtual ~PeptideStep() {}

    // The range given is the result from the previous step. Modifications to
    // the range for use by the next step should be stored back into the
    // provided range object. Similarly allow_detached indicates whether the
    // detached state is possible, and should be read and set appropriately.
    virtual void prune_forward(KDRange* range, bool* allow_detached) = 0;

    // The reverse of prune_forward. The range given is the result of calling
    // prune_backward() on the next step. Modifications to the range for use by
    // the previous step should be stored back into the range object. Similarly
    // allow_detached indicates whether the detached state is possible, and
    // should be read and set appropriately.
    virtual void prune_backward(KDRange* range, bool* allow_detached) = 0;
};

}  // namespace whatprot

#endif  // WHATPROT_HMM_STEP_PEPTIDE_STEP_H
