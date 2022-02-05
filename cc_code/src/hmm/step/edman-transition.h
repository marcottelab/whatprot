/******************************************************************************\
* Author: Matthew Beauregard Smith                                             *
* Affiliation: The University of Texas at Austin                               *
* Department: Oden Institute and Institute for Cellular and Molecular Biology  *
* PI: Edward Marcotte                                                          *
* Project: Protein Fluorosequencing                                            *
\******************************************************************************/

#ifndef WHATPROT_HMM_STEP_EDMAN_TRANSITION_H
#define WHATPROT_HMM_STEP_EDMAN_TRANSITION_H

// Local project headers:
#include "common/dye-track.h"
#include "hmm/state-vector/peptide-state-vector.h"
#include "hmm/step/peptide-step.h"
#include "parameterization/fit/sequencing-model-fitter.h"
#include "util/kd-range.h"

namespace whatprot {

class EdmanTransition : public PeptideStep {
public:
    EdmanTransition(double p_edman_failure,
                    const DyeSeq& dye_seq,
                    const DyeTrack& dye_track);
    EdmanTransition(const EdmanTransition& other);
    virtual ~EdmanTransition();
    void set_true_forward_range(const KDRange& range);
    void set_true_backward_range(const KDRange& range);
    virtual void prune_forward(KDRange* range, bool* allow_detached) override;
    virtual void prune_backward(KDRange* range, bool* allow_detached) override;
    virtual PeptideStateVector* forward(
            const PeptideStateVector& input,
            unsigned int* num_edmans) const override;
    virtual PeptideStateVector* backward(
            const PeptideStateVector& input,
            unsigned int* num_edmans) const override;
    virtual void improve_fit(const PeptideStateVector& forward_psv,
                             const PeptideStateVector& backward_psv,
                             const PeptideStateVector& next_backward_psv,
                             unsigned int num_edmans,
                             double probability,
                             SequencingModelFitter* fitter) const override;

    const DyeSeq& dye_seq;
    DyeTrack* dye_track;
    bool i_am_a_copy;
    double p_edman_failure;
    KDRange true_forward_range;  // shared with neighbors.
    KDRange safe_forward_range;  // used to make backward() easier.
    KDRange true_backward_range;  // shared with neighbors.
    KDRange safe_backward_range;  // used to make forward() easier.
};

}  // namespace whatprot

#endif  // WHATPROT_HMM_STEP_EDMAN_TRANSITION_H