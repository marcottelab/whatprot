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

namespace whatprot {

class EdmanTransition : public PeptideStep {
public:
    EdmanTransition(double p_edman_failure,
                    const DyeSeq& dye_seq,
                    const DyeTrack& dye_track);
    virtual void forward(unsigned int* num_edmans,
                         PeptideStateVector* psv) const override;
    virtual void backward(const PeptideStateVector& input,
                          unsigned int* num_edmans,
                          PeptideStateVector* output) const override;
    virtual void improve_fit(const PeptideStateVector& forward_psv,
                             const PeptideStateVector& backward_psv,
                             const PeptideStateVector& next_backward_psv,
                             unsigned int num_edmans,
                             double probability,
                             SequencingModelFitter* fitter) const override;

    DyeSeq dye_seq;
    DyeTrack dye_track;
    double p_edman_failure;
};

}  // namespace whatprot

#endif  // WHATPROT_HMM_STEP_EDMAN_TRANSITION_H