/******************************************************************************\
* Author: Matthew Beauregard Smith                                             *
* Affiliation: The University of Texas at Austin                               *
* Department: Oden Institute and Institute for Cellular and Molecular Biology  *
* PI: Edward Marcotte                                                          *
* Project: Protein Fluorosequencing                                            *
\******************************************************************************/

// Defining symbols from header:
#include "broken-n-transition.h"

// Local project headers:
#include "hmm/state-vector/peptide-state-vector.h"
#include "hmm/step/peptide-step.h"
#include "parameterization/fit/parameter-fitter.h"
#include "util/kd-range.h"

namespace whatprot {

BrokenNTransition::BrokenNTransition(double p_break_n) : p_break_n(p_break_n) {}

void BrokenNTransition::prune_forward(KDRange* range, bool* allow_detached) {
    pruned_range = *range;
}

void BrokenNTransition::prune_backward(KDRange* range, bool* allow_detached) {
    pruned_range = pruned_range.intersect(*range);
    *range = pruned_range;
}

PeptideStateVector* BrokenNTransition::forward(const PeptideStateVector& input,
                                               unsigned int* num_edmans) const {
    PeptideStateVector* output = new PeptideStateVector(pruned_range);
    ConstTensorIterator* tsr_in_itr = input.tensor.const_iterator(pruned_range);
    ConstTensorIterator* brkn_in_itr =
            input.broken_n_tensor.const_iterator(pruned_range);
    TensorIterator* tsr_out_itr = output->tensor.iterator(pruned_range);
    TensorIterator* brkn_out_itr =
            output->broken_n_tensor.iterator(pruned_range);
    while (!tsr_in_itr->done()) {
        *tsr_out_itr->get() = (1 - p_break_n) * (*tsr_in_itr->get());
        *brkn_out_itr->get() =
                (*brkn_in_itr->get()) + p_break_n * (*tsr_in_itr->get());
        tsr_in_itr->advance();
        brkn_in_itr->advance();
        tsr_out_itr->advance();
        brkn_out_itr->advance();
    }
    delete tsr_in_itr;
    delete brkn_in_itr;
    delete tsr_out_itr;
    delete brkn_out_itr;
    // Now we fix up the ranges, allow_detached, etc...
    output->range = pruned_range;
    output->allow_detached = input.allow_detached;
    if (output->allow_detached) {
        output->p_detached = input.p_detached;
    }
    return output;
}

PeptideStateVector* BrokenNTransition::backward(
        const PeptideStateVector& input, unsigned int* num_edmans) const {
    PeptideStateVector* output = new PeptideStateVector(pruned_range);
    ConstTensorIterator* tsr_in_itr = input.tensor.const_iterator(pruned_range);
    ConstTensorIterator* brkn_in_itr =
            input.broken_n_tensor.const_iterator(pruned_range);
    TensorIterator* tsr_out_itr = output->tensor.iterator(pruned_range);
    TensorIterator* brkn_out_itr =
            output->broken_n_tensor.iterator(pruned_range);
    while (!tsr_in_itr->done()) {
        *tsr_out_itr->get() = p_break_n * (*brkn_in_itr->get())
                              + (1 - p_break_n) * (*tsr_in_itr->get());
        *brkn_out_itr->get() = *brkn_in_itr->get();
        tsr_in_itr->advance();
        brkn_in_itr->advance();
        tsr_out_itr->advance();
        brkn_out_itr->advance();
    }
    delete tsr_in_itr;
    delete brkn_in_itr;
    delete tsr_out_itr;
    delete brkn_out_itr;
    // Now we fix up the ranges, allow_detached, etc...
    output->range = pruned_range;
    output->allow_detached = input.allow_detached;
    if (output->allow_detached) {
        output->p_detached = input.p_detached;
    }
    return output;
}

void BrokenNTransition::improve_fit(const PeptideStateVector& forward_psv,
                                    const PeptideStateVector& backward_psv,
                                    const PeptideStateVector& next_backward_psv,
                                    unsigned int num_edmans,
                                    double probability,
                                    ParameterFitter* fitter) const {
    ConstTensorIterator* f_itr =
            forward_psv.tensor.const_iterator(pruned_range);
    ConstTensorIterator* b_itr =
            backward_psv.tensor.const_iterator(pruned_range);
    ConstTensorIterator* nb_n_itr =
            next_backward_psv.broken_n_tensor.const_iterator(pruned_range);
    while (!f_itr->done()) {
        fitter->numerator +=
                (*f_itr->get()) * p_break_n * (*nb_n_itr->get()) / probability;
        fitter->denominator += (*f_itr->get()) * (*b_itr->get()) / probability;
        f_itr->advance();
        b_itr->advance();
        nb_n_itr->advance();
    }
    delete f_itr;
    delete b_itr;
    delete nb_n_itr;
}

}  // namespace whatprot
