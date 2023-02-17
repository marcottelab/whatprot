/******************************************************************************\
* Author: Matthew Beauregard Smith                                             *
* Affiliation: The University of Texas at Austin                               *
* Department: Oden Institute and Institute for Cellular and Molecular Biology  *
* PI: Edward Marcotte                                                          *
* Project: Protein Fluorosequencing                                            *
\******************************************************************************/

// Boost unit test framework (recommended to be the first include):
#include <boost/test/unit_test.hpp>

// File under test:
#include "block-transition.h"

// Local project headers:
#include "parameterization/fit/parameter-fitter.h"
#include "parameterization/fit/sequencing-model-fitter.h"
#include "util/kd-range.h"

namespace whatprot {

namespace {
using boost::unit_test::tolerance;
const double TOL = 0.000000001;
}  // namespace

// BrokenNTransition is abstract, so we need to override undefined methods.
class TestableBrokenNTransition : public BrokenNTransition {
public:
    TestableBrokenNTransition(double p_block)
            : BrokenNTransition(p_block) {}
    using BrokenNTransition::improve_fit;
    virtual void improve_fit(const PeptideStateVector& forward_psv,
                             const PeptideStateVector& backward_psv,
                             const PeptideStateVector& next_backward_psv,
                             unsigned int num_edmans,
                             double probability,
                             SequencingModelFitter* fitter) const override;
};

void TestableBrokenNTransition::improve_fit(
        const PeptideStateVector& forward_psv,
        const PeptideStateVector& backward_psv,
        const PeptideStateVector& next_backward_psv,
        unsigned int num_edmans,
        double probability,
        SequencingModelFitter* fitter) const {}

BOOST_AUTO_TEST_SUITE(hmm_suite)
BOOST_AUTO_TEST_SUITE(step_suite)
BOOST_AUTO_TEST_SUITE(broken_n_transition_suite)

BOOST_AUTO_TEST_CASE(forward_test) {
    double p_block = 0.07;
    TestableBrokenNTransition bnt(p_block);
    bnt.pruned_range.min = {0, 0};
    bnt.pruned_range.max = {1, 2};
    unsigned int order = 2;
    unsigned int* shape = new unsigned int[order];
    shape[0] = 1;
    shape[1] = 2;
    PeptideStateVector psv1(order, shape);
    delete[] shape;
    psv1.tensor[{0, 0}] = 0.3;
    psv1.tensor[{0, 1}] = 0.7;
    psv1.broken_n_tensor[{0, 0}] = 0.03;
    psv1.broken_n_tensor[{0, 1}] = 0.07;
    psv1.allow_detached = true;
    psv1.p_detached = 0.123;
    unsigned int edmans = 0;
    PeptideStateVector* psv2 = bnt.forward(psv1, &edmans);
    BOOST_TEST((psv2->tensor[{0, 0}]) == (1 - p_block) * 0.3);
    BOOST_TEST((psv2->tensor[{0, 1}]) == (1 - p_block) * 0.7);
    BOOST_TEST((psv2->broken_n_tensor[{0, 0}]) == 0.03 + p_block * 0.3);
    BOOST_TEST((psv2->broken_n_tensor[{0, 1}]) == 0.07 + p_block * 0.7);
    BOOST_TEST(psv2->range.min[0] == 0u);
    BOOST_TEST(psv2->range.min[1] == 0u);
    BOOST_TEST(psv2->range.max[0] == 1u);
    BOOST_TEST(psv2->range.max[1] == 2u);
    BOOST_TEST(psv2->allow_detached == true);
    BOOST_TEST(psv2->p_detached == 0.123);
    delete psv2;
}

BOOST_AUTO_TEST_CASE(backward_test) {
    double p_block = 0.07;
    TestableBrokenNTransition bnt(p_block);
    bnt.pruned_range.min = {0, 0};
    bnt.pruned_range.max = {1, 2};
    unsigned int order = 2;
    unsigned int* shape = new unsigned int[order];
    shape[0] = 1;
    shape[1] = 2;
    PeptideStateVector psv1(order, shape);
    delete[] shape;
    psv1.tensor[{0, 0}] = 0.3;
    psv1.tensor[{0, 1}] = 0.7;
    psv1.broken_n_tensor[{0, 0}] = 0.03;
    psv1.broken_n_tensor[{0, 1}] = 0.07;
    psv1.allow_detached = true;
    psv1.p_detached = 0.123;
    unsigned int edmans = 0;
    PeptideStateVector* psv2 = bnt.backward(psv1, &edmans);
    BOOST_TEST((psv2->tensor[{0, 0}])
               == p_block * 0.03 + (1 - p_block) * 0.3);
    BOOST_TEST((psv2->tensor[{0, 1}])
               == p_block * 0.07 + (1 - p_block) * 0.7);
    BOOST_TEST((psv2->broken_n_tensor[{0, 0}]) == 0.03);
    BOOST_TEST((psv2->broken_n_tensor[{0, 1}]) == 0.07);
    BOOST_TEST(psv2->range.min[0] == 0u);
    BOOST_TEST(psv2->range.min[1] == 0u);
    BOOST_TEST(psv2->range.max[0] == 1u);
    BOOST_TEST(psv2->range.max[1] == 2u);
    BOOST_TEST(psv2->allow_detached == true);
    BOOST_TEST(psv2->p_detached == 0.123);
    delete psv2;
}

BOOST_AUTO_TEST_CASE(improve_fit_test) {
    double p_block = 0.07;
    TestableBrokenNTransition bnt(p_block);
    bnt.pruned_range.min = {0, 0};
    bnt.pruned_range.max = {1, 2};
    unsigned int order = 2;
    unsigned int* shape = new unsigned int[order];
    shape[0] = 1;
    shape[1] = 2;
    PeptideStateVector fpsv(order, shape);
    fpsv.tensor[{0, 0}] = 0.31;
    fpsv.tensor[{0, 1}] = 0.71;
    fpsv.broken_n_tensor[{0, 0}] = 0.031;
    fpsv.broken_n_tensor[{0, 1}] = 0.071;
    PeptideStateVector bpsv(order, shape);
    bpsv.tensor[{0, 0}] = 0.32;
    bpsv.tensor[{0, 1}] = 0.72;
    bpsv.broken_n_tensor[{0, 0}] = 0.032;
    bpsv.broken_n_tensor[{0, 1}] = 0.072;
    PeptideStateVector nbpsv(order, shape);
    nbpsv.tensor[{0, 0}] = 0.33;
    nbpsv.tensor[{0, 1}] = 0.73;
    nbpsv.broken_n_tensor[{0, 0}] = 0.033;
    nbpsv.broken_n_tensor[{0, 1}] = 0.073;
    delete[] shape;
    unsigned int edmans = 0;
    double probability = 1.0;
    ParameterFitter pf;
    bnt.improve_fit(fpsv, bpsv, nbpsv, edmans, probability, &pf);
    BOOST_TEST(pf.get()
               == (0.31 * p_block * 0.033 + 0.71 * p_block * 0.073)
                          / (0.31 * 0.32 + 0.71 * 0.72));
}

BOOST_AUTO_TEST_SUITE_END()  // broken_n_transition_suite
BOOST_AUTO_TEST_SUITE_END()  // step_suite
BOOST_AUTO_TEST_SUITE_END()  // hmm_suite

}  // namespace whatprot
