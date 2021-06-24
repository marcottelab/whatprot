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
#include "binomial-transition.h"

// Local project headers:
#include "parameterization/fit/parameter-fitter.h"
#include "parameterization/fit/sequencing-model-fitter.h"

namespace whatprot {

namespace {
using boost::unit_test::tolerance;
const double TOL = 0.000000001;
}  // namespace

// BinomialTransition is abstract, so we need to override undefined methods.
class TestableBinomialTransition : public BinomialTransition {
public:
    TestableBinomialTransition(double q, int channel)
            : BinomialTransition(q, channel) {}
    using BinomialTransition::improve_fit;
    virtual void improve_fit(const PeptideStateVector& forward_psv,
                             const PeptideStateVector& backward_psv,
                             const PeptideStateVector& next_backward_psv,
                             unsigned int num_edmans,
                             double probability,
                             SequencingModelFitter* fitter) const override;
};

void TestableBinomialTransition::improve_fit(
        const PeptideStateVector& forward_psv,
        const PeptideStateVector& backward_psv,
        const PeptideStateVector& next_backward_psv,
        unsigned int num_edmans,
        double probability,
        SequencingModelFitter* fitter) const {}

BOOST_AUTO_TEST_SUITE(hmm_suite)
BOOST_AUTO_TEST_SUITE(step_suite)
BOOST_AUTO_TEST_SUITE(binomial_transition_suite)

BOOST_AUTO_TEST_CASE(constructor_test, *tolerance(TOL)) {
    double q = 0.2;
    int channel = -1;  // can be ignored for this test.
    TestableBinomialTransition bt(q, channel);
    BOOST_TEST(bt.q == q);
    BOOST_TEST(bt.channel == channel);
}

BOOST_AUTO_TEST_CASE(reserve_zero_test, *tolerance(TOL)) {
    double q = 0.2;
    int channel = -1;  // can be ignored for this test.
    TestableBinomialTransition bt(q, channel);
    bt.reserve(0);
    BOOST_TEST(bt.prob(0, 0) == 1.0);
}

BOOST_AUTO_TEST_CASE(reserve_one_test, *tolerance(TOL)) {
    double q = 0.2;
    double p = 0.8;
    int channel = -1;  // can be ignored for this test.
    TestableBinomialTransition bt(q, channel);
    bt.reserve(1);
    BOOST_TEST(bt.prob(0, 0) == 1.0);
    BOOST_TEST(bt.prob(1, 0) == q);
    BOOST_TEST(bt.prob(1, 1) == p);
}

BOOST_AUTO_TEST_CASE(reserve_two_test, *tolerance(TOL)) {
    double q = 0.2;
    double p = 0.8;
    int channel = -1;  // can be ignored for this test.
    TestableBinomialTransition bt(q, channel);
    bt.reserve(2);
    BOOST_TEST(bt.prob(0, 0) == 1.0);
    BOOST_TEST(bt.prob(1, 0) == q);
    BOOST_TEST(bt.prob(1, 1) == p);
    BOOST_TEST(bt.prob(2, 0) == q * q);
    BOOST_TEST(bt.prob(2, 1) == 2 * q * p);
    BOOST_TEST(bt.prob(2, 2) == p * p);
}

BOOST_AUTO_TEST_CASE(reserve_three_test, *tolerance(TOL)) {
    double q = 0.2;
    double p = 0.8;
    int channel = -1;  // can be ignored for this test.
    TestableBinomialTransition bt(q, channel);
    bt.reserve(3);
    BOOST_TEST(bt.prob(0, 0) == 1.0);
    BOOST_TEST(bt.prob(1, 0) == q);
    BOOST_TEST(bt.prob(1, 1) == p);
    BOOST_TEST(bt.prob(2, 0) == q * q);
    BOOST_TEST(bt.prob(2, 1) == 2 * q * p);
    BOOST_TEST(bt.prob(2, 2) == p * p);
    BOOST_TEST(bt.prob(3, 0) == q * q * q);
    BOOST_TEST(bt.prob(3, 1) == 3 * q * q * p);
    BOOST_TEST(bt.prob(3, 2) == 3 * q * p * p);
    BOOST_TEST(bt.prob(3, 3) == p * p * p);
}

BOOST_AUTO_TEST_CASE(reserve_no_shrink_test, *tolerance(TOL)) {
    double q = 0.2;
    double p = 0.8;
    int channel = -1;  // can be ignored for this test.
    TestableBinomialTransition bt(q, channel);
    bt.reserve(3);
    bt.reserve(1);  // This operation should be silently ignored.
    BOOST_TEST(bt.prob(0, 0) == 1.0);
    BOOST_TEST(bt.prob(1, 0) == q);
    BOOST_TEST(bt.prob(1, 1) == p);
    BOOST_TEST(bt.prob(2, 0) == q * q);
    BOOST_TEST(bt.prob(2, 1) == 2 * q * p);
    BOOST_TEST(bt.prob(2, 2) == p * p);
    BOOST_TEST(bt.prob(3, 0) == q * q * q);
    BOOST_TEST(bt.prob(3, 1) == 3 * q * q * p);
    BOOST_TEST(bt.prob(3, 2) == 3 * q * p * p);
    BOOST_TEST(bt.prob(3, 3) == p * p * p);
}

BOOST_AUTO_TEST_CASE(forward_in_place_trivial_test, *tolerance(TOL)) {
    double q = 0.05;
    int channel = 0;
    TestableBinomialTransition bt(q, channel);
    bt.reserve(0);
    unsigned int order = 2;
    unsigned int* shape = new unsigned int[order];
    shape[0] = 1;
    shape[1] = 1;
    PeptideStateVector psv(order, shape);
    delete[] shape;
    psv.tensor[{0, 0}] = 1.0;
    unsigned int edmans = 0;
    bt.forward(&edmans, &psv);
    BOOST_TEST((psv.tensor[{0, 0}]) == 1.0);
}

BOOST_AUTO_TEST_CASE(forward_in_place_basic_transition_test, *tolerance(TOL)) {
    double q = 0.05;
    double p = 0.95;
    int channel = 0;
    TestableBinomialTransition bt(q, channel);
    bt.reserve(1);
    unsigned int order = 2;
    unsigned int* shape = new unsigned int[order];
    shape[0] = 1;
    shape[1] = 2;
    PeptideStateVector psv(order, shape);
    delete[] shape;
    psv.tensor[{0, 0}] = 0.3;
    psv.tensor[{0, 1}] = 0.7;
    unsigned int edmans = 0;
    bt.forward(&edmans, &psv);
    BOOST_TEST((psv.tensor[{0, 0}]) == 0.3 + 0.7 * q);
    BOOST_TEST((psv.tensor[{0, 1}]) == 0.7 * p);
}

BOOST_AUTO_TEST_CASE(forward_in_place_bigger_transition_test, *tolerance(TOL)) {
    double q = 0.05;
    double p = 0.95;
    int channel = 0;
    TestableBinomialTransition bt(q, channel);
    bt.reserve(2);
    unsigned int order = 2;
    unsigned int* shape = new unsigned int[order];
    shape[0] = 1;
    shape[1] = 3;
    PeptideStateVector psv(order, shape);
    delete[] shape;
    psv.tensor[{0, 0}] = 0.2;
    psv.tensor[{0, 1}] = 0.3;
    psv.tensor[{0, 2}] = 0.7;
    unsigned int edmans = 0;
    bt.forward(&edmans, &psv);
    BOOST_TEST((psv.tensor[{0, 0}]) == 0.2 + 0.3 * q + 0.7 * q * q);
    BOOST_TEST((psv.tensor[{0, 1}]) == 0.3 * p + 0.7 * 2 * q * p);
    BOOST_TEST((psv.tensor[{0, 2}]) == 0.7 * p * p);
}

BOOST_AUTO_TEST_CASE(forward_in_place_multiple_edmans_test, *tolerance(TOL)) {
    double q = 0.05;
    double p = 0.95;
    int channel = 0;
    TestableBinomialTransition bt(q, channel);
    bt.reserve(2);
    unsigned int order = 2;
    unsigned int* shape = new unsigned int[order];
    shape[0] = 3;
    shape[1] = 2;
    PeptideStateVector psv(order, shape);
    delete[] shape;
    psv.tensor[{0, 0}] = 0.2;
    psv.tensor[{0, 1}] = 0.8;
    psv.tensor[{1, 0}] = 0.3;
    psv.tensor[{1, 1}] = 0.7;
    psv.tensor[{2, 0}] = 0.4;
    psv.tensor[{2, 1}] = 0.6;
    unsigned int edmans = 2;
    bt.forward(&edmans, &psv);
    BOOST_TEST((psv.tensor[{0, 0}]) == 0.2 + 0.8 * q);
    BOOST_TEST((psv.tensor[{0, 1}]) == 0.8 * p);
    BOOST_TEST((psv.tensor[{1, 0}]) == 0.3 + 0.7 * q);
    BOOST_TEST((psv.tensor[{1, 1}]) == 0.7 * p);
    BOOST_TEST((psv.tensor[{2, 0}]) == 0.4 + 0.6 * q);
    BOOST_TEST((psv.tensor[{2, 1}]) == 0.6 * p);
}

BOOST_AUTO_TEST_CASE(forward_in_place_other_dye_colors_test, *tolerance(TOL)) {
    double q = 0.05;
    double p = 0.95;
    int channel = 1;  // corresponds to 2nd dim of tensor
    TestableBinomialTransition bt(q, channel);
    bt.reserve(1);
    unsigned int order = 4;
    unsigned int* shape = new unsigned int[order];
    shape[0] = 1;
    shape[1] = 2;
    shape[2] = 2;
    shape[3] = 2;
    PeptideStateVector psv(order, shape);
    delete[] shape;
    psv.tensor[{0, 0, 0, 0}] = 0.1;
    psv.tensor[{0, 0, 0, 1}] = 0.2;
    psv.tensor[{0, 0, 1, 0}] = 0.3;
    psv.tensor[{0, 0, 1, 1}] = 0.4;
    psv.tensor[{0, 1, 0, 0}] = 0.5;
    psv.tensor[{0, 1, 0, 1}] = 0.6;
    psv.tensor[{0, 1, 1, 0}] = 0.7;
    psv.tensor[{0, 1, 1, 1}] = 0.8;
    unsigned int edmans = 0;
    bt.forward(&edmans, &psv);
    BOOST_TEST((psv.tensor[{0, 0, 0, 0}]) == 0.1 + 0.3 * q);
    BOOST_TEST((psv.tensor[{0, 0, 0, 1}]) == 0.2 + 0.4 * q);
    BOOST_TEST((psv.tensor[{0, 0, 1, 0}]) == 0.3 * p);
    BOOST_TEST((psv.tensor[{0, 0, 1, 1}]) == 0.4 * p);
    BOOST_TEST((psv.tensor[{0, 1, 0, 0}]) == 0.5 + 0.7 * q);
    BOOST_TEST((psv.tensor[{0, 1, 0, 1}]) == 0.6 + 0.8 * q);
    BOOST_TEST((psv.tensor[{0, 1, 1, 0}]) == 0.7 * p);
    BOOST_TEST((psv.tensor[{0, 1, 1, 1}]) == 0.8 * p);
}

BOOST_AUTO_TEST_CASE(backward_in_place_trivial_test, *tolerance(TOL)) {
    double q = 0.05;
    int channel = 0;
    TestableBinomialTransition bt(q, channel);
    bt.reserve(0);
    unsigned int order = 2;
    unsigned int* shape = new unsigned int[order];
    shape[0] = 1;
    shape[1] = 1;
    PeptideStateVector psv(order, shape);
    delete[] shape;
    psv.tensor[{0, 0}] = 1.0;
    unsigned int edmans = 0;
    bt.backward(psv, &edmans, &psv);
    BOOST_TEST((psv.tensor[{0, 0}]) == 1.0);
}

BOOST_AUTO_TEST_CASE(backward_new_tsr_trivial_test, *tolerance(TOL)) {
    double q = 0.05;
    int channel = 0;
    TestableBinomialTransition bt(q, channel);
    bt.reserve(0);
    unsigned int order = 2;
    unsigned int* shape = new unsigned int[order];
    shape[0] = 1;
    shape[1] = 1;
    PeptideStateVector psv1(order, shape);
    PeptideStateVector psv2(order, shape);
    delete[] shape;
    psv1.tensor[{0, 0}] = 1.0;
    unsigned int edmans = 0;
    bt.backward(psv1, &edmans, &psv2);
    BOOST_TEST((psv2.tensor[{0, 0}]) == 1.0);
}

BOOST_AUTO_TEST_CASE(backward_in_place_basic_transition_test, *tolerance(TOL)) {
    double q = 0.05;
    double p = 0.95;
    int channel = 0;
    TestableBinomialTransition bt(q, channel);
    bt.reserve(1);
    unsigned int order = 2;
    unsigned int* shape = new unsigned int[order];
    shape[0] = 1;
    shape[1] = 2;
    PeptideStateVector psv(order, shape);
    delete[] shape;
    psv.tensor[{0, 0}] = 0.3;
    psv.tensor[{0, 1}] = 0.7;
    unsigned int edmans = 0;
    bt.backward(psv, &edmans, &psv);
    BOOST_TEST((psv.tensor[{0, 0}]) == 0.3);
    BOOST_TEST((psv.tensor[{0, 1}]) == q * 0.3 + p * 0.7);
}

BOOST_AUTO_TEST_CASE(backward_new_tsr_basic_transition_test, *tolerance(TOL)) {
    double q = 0.05;
    double p = 0.95;
    int channel = 0;
    TestableBinomialTransition bt(q, channel);
    bt.reserve(1);
    unsigned int order = 2;
    unsigned int* shape = new unsigned int[order];
    shape[0] = 1;
    shape[1] = 2;
    PeptideStateVector psv1(order, shape);
    PeptideStateVector psv2(order, shape);
    delete[] shape;
    psv1.tensor[{0, 0}] = 0.3;
    psv1.tensor[{0, 1}] = 0.7;
    unsigned int edmans = 0;
    bt.backward(psv1, &edmans, &psv2);
    BOOST_TEST((psv2.tensor[{0, 0}]) == 0.3);
    BOOST_TEST((psv2.tensor[{0, 1}]) == q * 0.3 + p * 0.7);
}

BOOST_AUTO_TEST_CASE(backward_in_place_bigger_transition_test,
                     *tolerance(TOL)) {
    double q = 0.05;
    double p = 0.95;
    int channel = 0;
    TestableBinomialTransition bt(q, channel);
    bt.reserve(2);
    unsigned int order = 2;
    unsigned int* shape = new unsigned int[order];
    shape[0] = 1;
    shape[1] = 3;
    PeptideStateVector psv(order, shape);
    delete[] shape;
    psv.tensor[{0, 0}] = 0.2;
    psv.tensor[{0, 1}] = 0.3;
    psv.tensor[{0, 2}] = 0.7;
    unsigned int edmans = 0;
    bt.backward(psv, &edmans, &psv);
    BOOST_TEST((psv.tensor[{0, 0}]) == 0.2);
    BOOST_TEST((psv.tensor[{0, 1}]) == q * 0.2 + p * 0.3);
    BOOST_TEST((psv.tensor[{0, 2}])
               == q * q * 0.2 + 2 * q * p * 0.3 + p * p * 0.7);
}

BOOST_AUTO_TEST_CASE(backward_new_tsr_bigger_transition_test, *tolerance(TOL)) {
    double q = 0.05;
    double p = 0.95;
    int channel = 0;
    TestableBinomialTransition bt(q, channel);
    bt.reserve(2);
    unsigned int order = 2;
    unsigned int* shape = new unsigned int[order];
    shape[0] = 1;
    shape[1] = 3;
    PeptideStateVector psv1(order, shape);
    PeptideStateVector psv2(order, shape);
    delete[] shape;
    psv1.tensor[{0, 0}] = 0.2;
    psv1.tensor[{0, 1}] = 0.3;
    psv1.tensor[{0, 2}] = 0.7;
    unsigned int edmans = 0;
    bt.backward(psv1, &edmans, &psv2);
    BOOST_TEST((psv2.tensor[{0, 0}]) == 0.2);
    BOOST_TEST((psv2.tensor[{0, 1}]) == q * 0.2 + p * 0.3);
    BOOST_TEST((psv2.tensor[{0, 2}])
               == q * q * 0.2 + 2 * q * p * 0.3 + p * p * 0.7);
}

BOOST_AUTO_TEST_CASE(backward_in_place_multiple_edmans_test, *tolerance(TOL)) {
    double q = 0.05;
    double p = 0.95;
    int channel = 0;
    TestableBinomialTransition bt(q, channel);
    bt.reserve(2);
    unsigned int order = 2;
    unsigned int* shape = new unsigned int[order];
    shape[0] = 3;
    shape[1] = 2;
    PeptideStateVector psv(order, shape);
    delete[] shape;
    psv.tensor[{0, 0}] = 0.2;
    psv.tensor[{0, 1}] = 0.8;
    psv.tensor[{1, 0}] = 0.3;
    psv.tensor[{1, 1}] = 0.7;
    psv.tensor[{2, 0}] = 0.4;
    psv.tensor[{2, 1}] = 0.6;
    unsigned int edmans = 2;
    bt.backward(psv, &edmans, &psv);
    BOOST_TEST((psv.tensor[{0, 0}]) == 0.2);
    BOOST_TEST((psv.tensor[{0, 1}]) == q * 0.2 + p * 0.8);
    BOOST_TEST((psv.tensor[{1, 0}]) == 0.3);
    BOOST_TEST((psv.tensor[{1, 1}]) == q * 0.3 + p * 0.7);
    BOOST_TEST((psv.tensor[{2, 0}]) == 0.4);
    BOOST_TEST((psv.tensor[{2, 1}]) == q * 0.4 + p * 0.6);
}

BOOST_AUTO_TEST_CASE(backward_new_tsr_multiple_edmans_test, *tolerance(TOL)) {
    double q = 0.05;
    double p = 0.95;
    int channel = 0;
    TestableBinomialTransition bt(q, channel);
    bt.reserve(2);
    unsigned int order = 2;
    unsigned int* shape = new unsigned int[order];
    shape[0] = 3;
    shape[1] = 2;
    PeptideStateVector psv1(order, shape);
    PeptideStateVector psv2(order, shape);
    delete[] shape;
    psv1.tensor[{0, 0}] = 0.2;
    psv1.tensor[{0, 1}] = 0.8;
    psv1.tensor[{1, 0}] = 0.3;
    psv1.tensor[{1, 1}] = 0.7;
    psv1.tensor[{2, 0}] = 0.4;
    psv1.tensor[{2, 1}] = 0.6;
    unsigned int edmans = 2;
    bt.backward(psv1, &edmans, &psv2);
    BOOST_TEST((psv2.tensor[{0, 0}]) == 0.2);
    BOOST_TEST((psv2.tensor[{0, 1}]) == q * 0.2 + p * 0.8);
    BOOST_TEST((psv2.tensor[{1, 0}]) == 0.3);
    BOOST_TEST((psv2.tensor[{1, 1}]) == q * 0.3 + p * 0.7);
    BOOST_TEST((psv2.tensor[{2, 0}]) == 0.4);
    BOOST_TEST((psv2.tensor[{2, 1}]) == q * 0.4 + p * 0.6);
}

BOOST_AUTO_TEST_CASE(backward_in_place_other_dye_colors_test, *tolerance(TOL)) {
    double q = 0.05;
    double p = 0.95;
    int channel = 1;  // corresponds to 2nd dim of tensor
    TestableBinomialTransition bt(q, channel);
    bt.reserve(1);
    unsigned int order = 4;
    unsigned int* shape = new unsigned int[order];
    shape[0] = 1;
    shape[1] = 2;
    shape[2] = 2;
    shape[3] = 2;
    PeptideStateVector psv(order, shape);
    delete[] shape;
    psv.tensor[{0, 0, 0, 0}] = 0.1;
    psv.tensor[{0, 0, 0, 1}] = 0.2;
    psv.tensor[{0, 0, 1, 0}] = 0.3;
    psv.tensor[{0, 0, 1, 1}] = 0.4;
    psv.tensor[{0, 1, 0, 0}] = 0.5;
    psv.tensor[{0, 1, 0, 1}] = 0.6;
    psv.tensor[{0, 1, 1, 0}] = 0.7;
    psv.tensor[{0, 1, 1, 1}] = 0.8;
    unsigned int edmans = 0;
    bt.backward(psv, &edmans, &psv);
    BOOST_TEST((psv.tensor[{0, 0, 0, 0}]) == 0.1);
    BOOST_TEST((psv.tensor[{0, 0, 0, 1}]) == 0.2);
    BOOST_TEST((psv.tensor[{0, 0, 1, 0}]) == q * 0.1 + p * 0.3);
    BOOST_TEST((psv.tensor[{0, 0, 1, 1}]) == q * 0.2 + p * 0.4);
    BOOST_TEST((psv.tensor[{0, 1, 0, 0}]) == 0.5);
    BOOST_TEST((psv.tensor[{0, 1, 0, 1}]) == 0.6);
    BOOST_TEST((psv.tensor[{0, 1, 1, 0}]) == q * 0.5 + p * 0.7);
    BOOST_TEST((psv.tensor[{0, 1, 1, 1}]) == q * 0.6 + p * 0.8);
}

BOOST_AUTO_TEST_CASE(backward_new_tsr_other_dye_colors_test, *tolerance(TOL)) {
    double q = 0.05;
    double p = 0.95;
    int channel = 1;  // corresponds to 2nd dim of tensor
    TestableBinomialTransition bt(q, channel);
    bt.reserve(1);
    unsigned int order = 4;
    unsigned int* shape = new unsigned int[order];
    shape[0] = 1;
    shape[1] = 2;
    shape[2] = 2;
    shape[3] = 2;
    PeptideStateVector psv1(order, shape);
    PeptideStateVector psv2(order, shape);
    delete[] shape;
    psv1.tensor[{0, 0, 0, 0}] = 0.1;
    psv1.tensor[{0, 0, 0, 1}] = 0.2;
    psv1.tensor[{0, 0, 1, 0}] = 0.3;
    psv1.tensor[{0, 0, 1, 1}] = 0.4;
    psv1.tensor[{0, 1, 0, 0}] = 0.5;
    psv1.tensor[{0, 1, 0, 1}] = 0.6;
    psv1.tensor[{0, 1, 1, 0}] = 0.7;
    psv1.tensor[{0, 1, 1, 1}] = 0.8;
    unsigned int edmans = 0;
    bt.backward(psv1, &edmans, &psv2);
    BOOST_TEST((psv2.tensor[{0, 0, 0, 0}]) == 0.1);
    BOOST_TEST((psv2.tensor[{0, 0, 0, 1}]) == 0.2);
    BOOST_TEST((psv2.tensor[{0, 0, 1, 0}]) == q * 0.1 + p * 0.3);
    BOOST_TEST((psv2.tensor[{0, 0, 1, 1}]) == q * 0.2 + p * 0.4);
    BOOST_TEST((psv2.tensor[{0, 1, 0, 0}]) == 0.5);
    BOOST_TEST((psv2.tensor[{0, 1, 0, 1}]) == 0.6);
    BOOST_TEST((psv2.tensor[{0, 1, 1, 0}]) == q * 0.5 + p * 0.7);
    BOOST_TEST((psv2.tensor[{0, 1, 1, 1}]) == q * 0.6 + p * 0.8);
}

BOOST_AUTO_TEST_CASE(improve_fit_trivial_test, *tolerance(TOL)) {
    double q = 0.05;
    int channel = 0;
    TestableBinomialTransition bt(q, channel);
    bt.reserve(0);
    unsigned int order = 2;
    unsigned int* shape = new unsigned int[order];
    shape[0] = 1;
    shape[1] = 1;
    PeptideStateVector fpsv(order, shape);
    fpsv.tensor[{0, 0}] = 1.0;
    PeptideStateVector bpsv(order, shape);
    bpsv.tensor[{0, 0}] = 1.0;
    PeptideStateVector nbpsv(order, shape);
    nbpsv.tensor[{0, 0}] = 1.0;
    delete[] shape;
    unsigned int edmans = 0;
    double probability = 1.0;
    ParameterFitter pf;
    bt.improve_fit(fpsv, bpsv, nbpsv, edmans, probability, &pf);
    BOOST_TEST(pf.numerator == 0.0);
    BOOST_TEST(pf.denominator == 0.0);
}

BOOST_AUTO_TEST_CASE(improve_fit_basic_test, *tolerance(TOL)) {
    double q = 0.05;
    int channel = 0;
    TestableBinomialTransition bt(q, channel);
    bt.reserve(1);
    unsigned int order = 2;
    unsigned int* shape = new unsigned int[order];
    shape[0] = 1;
    shape[1] = 2;
    PeptideStateVector fpsv(order, shape);
    fpsv.tensor[{0, 0}] = 0.31;
    fpsv.tensor[{0, 1}] = 0.71;
    PeptideStateVector bpsv(order, shape);
    bpsv.tensor[{0, 0}] = 0.32;
    bpsv.tensor[{0, 1}] = 0.72;
    PeptideStateVector nbpsv(order, shape);
    nbpsv.tensor[{0, 0}] = 0.33;
    nbpsv.tensor[{0, 1}] = 0.73;
    delete[] shape;
    unsigned int edmans = 0;
    double probability = 1.0;
    ParameterFitter pf;
    bt.improve_fit(fpsv, bpsv, nbpsv, edmans, probability, &pf);
    BOOST_TEST(pf.get() == (0.71 * q * 0.33) / (0.71 * 0.72));
}

BOOST_AUTO_TEST_CASE(improve_fit_bigger_test, *tolerance(TOL)) {
    double q = 0.05;
    double p = 0.95;
    int channel = 0;
    TestableBinomialTransition bt(q, channel);
    bt.reserve(2);
    unsigned int order = 2;
    unsigned int* shape = new unsigned int[order];
    shape[0] = 1;
    shape[1] = 3;
    PeptideStateVector fpsv(order, shape);
    fpsv.tensor[{0, 0}] = 0.31;
    fpsv.tensor[{0, 1}] = 0.71;
    fpsv.tensor[{0, 2}] = 0.91;
    PeptideStateVector bpsv(order, shape);
    bpsv.tensor[{0, 0}] = 0.32;
    bpsv.tensor[{0, 1}] = 0.72;
    bpsv.tensor[{0, 2}] = 0.92;
    PeptideStateVector nbpsv(order, shape);
    nbpsv.tensor[{0, 0}] = 0.33;
    nbpsv.tensor[{0, 1}] = 0.73;
    nbpsv.tensor[{0, 2}] = 0.93;
    delete[] shape;
    unsigned int edmans = 0;
    double probability = 1.0;
    ParameterFitter pf;
    bt.improve_fit(fpsv, bpsv, nbpsv, edmans, probability, &pf);
    BOOST_TEST(pf.get()
               == (0.71 * q * 0.33 + 0.91 * (q * p * 2.0) * 0.73
                   + 0.91 * (q * q) * 0.33 * 2.0)
                          / (0.71 * 0.72 + 0.91 * 0.92 * 2.0));
}

BOOST_AUTO_TEST_CASE(improve_fit_multiple_edmans_test, *tolerance(TOL)) {
    double q = 0.05;
    int channel = 0;
    TestableBinomialTransition bt(q, channel);
    bt.reserve(1);
    unsigned int order = 2;
    unsigned int* shape = new unsigned int[order];
    shape[0] = 2;
    shape[1] = 2;
    PeptideStateVector fpsv(order, shape);
    fpsv.tensor[{0, 0}] = 0.31;
    fpsv.tensor[{0, 1}] = 0.71;
    fpsv.tensor[{1, 0}] = 0.41;
    fpsv.tensor[{1, 1}] = 0.81;
    PeptideStateVector bpsv(order, shape);
    bpsv.tensor[{0, 0}] = 0.32;
    bpsv.tensor[{0, 1}] = 0.72;
    bpsv.tensor[{1, 0}] = 0.42;
    bpsv.tensor[{1, 1}] = 0.82;
    PeptideStateVector nbpsv(order, shape);
    nbpsv.tensor[{0, 0}] = 0.33;
    nbpsv.tensor[{0, 1}] = 0.73;
    nbpsv.tensor[{1, 0}] = 0.43;
    nbpsv.tensor[{1, 1}] = 0.83;
    delete[] shape;
    unsigned int edmans = 1;
    double probability = 1.0;
    ParameterFitter pf;
    bt.improve_fit(fpsv, bpsv, nbpsv, edmans, probability, &pf);
    BOOST_TEST(pf.get()
               == (0.71 * q * 0.33 + 0.81 * q * 0.43)
                          / (0.71 * 0.72 + 0.81 * 0.82));
}

BOOST_AUTO_TEST_CASE(improve_fit_other_dye_color_test, *tolerance(TOL)) {
    double q = 0.05;
    int channel = 0;
    TestableBinomialTransition bt(q, channel);
    bt.reserve(1);
    unsigned int order = 3;
    unsigned int* shape = new unsigned int[order];
    shape[0] = 1;
    shape[1] = 2;
    shape[2] = 2;
    PeptideStateVector fpsv(order, shape);
    fpsv.tensor[{0, 0, 0}] = 0.31;
    fpsv.tensor[{0, 0, 1}] = 0.231;
    fpsv.tensor[{0, 1, 0}] = 0.71;
    fpsv.tensor[{0, 1, 1}] = 0.271;
    PeptideStateVector bpsv(order, shape);
    bpsv.tensor[{0, 0, 0}] = 0.32;
    bpsv.tensor[{0, 0, 1}] = 0.232;
    bpsv.tensor[{0, 1, 0}] = 0.72;
    bpsv.tensor[{0, 1, 1}] = 0.272;
    PeptideStateVector nbpsv(order, shape);
    nbpsv.tensor[{0, 0, 0}] = 0.33;
    nbpsv.tensor[{0, 0, 1}] = 0.233;
    nbpsv.tensor[{0, 1, 0}] = 0.73;
    nbpsv.tensor[{0, 1, 1}] = 0.273;
    delete[] shape;
    unsigned int edmans = 0;
    double probability = 1.0;
    ParameterFitter pf;
    bt.improve_fit(fpsv, bpsv, nbpsv, edmans, probability, &pf);
    BOOST_TEST(pf.get()
               == (0.71 * q * 0.33 + 0.271 * q * 0.233)
                          / (0.71 * 0.72 + 0.271 * 0.272));
}

BOOST_AUTO_TEST_CASE(improve_fit_different_probability_test, *tolerance(TOL)) {
    double q = 0.05;
    int channel = 0;
    TestableBinomialTransition bt(q, channel);
    bt.reserve(1);
    unsigned int order = 2;
    unsigned int* shape = new unsigned int[order];
    shape[0] = 1;
    shape[1] = 2;
    PeptideStateVector fpsv(order, shape);
    fpsv.tensor[{0, 0}] = 0.31;
    fpsv.tensor[{0, 1}] = 0.71;
    PeptideStateVector bpsv(order, shape);
    bpsv.tensor[{0, 0}] = 0.32;
    bpsv.tensor[{0, 1}] = 0.72;
    PeptideStateVector nbpsv(order, shape);
    nbpsv.tensor[{0, 0}] = 0.33;
    nbpsv.tensor[{0, 1}] = 0.73;
    delete[] shape;
    unsigned int edmans = 0;
    double probability = 0.123456789;
    ParameterFitter pf;
    bt.improve_fit(fpsv, bpsv, nbpsv, edmans, probability, &pf);
    BOOST_TEST(pf.get() == (0.71 * q * 0.33) / (0.71 * 0.72));
}

BOOST_AUTO_TEST_CASE(improve_fit_twice_test, *tolerance(TOL)) {
    double q = 0.05;
    int channel = 0;
    TestableBinomialTransition bt(q, channel);
    bt.reserve(1);
    unsigned int order = 2;
    unsigned int* shape = new unsigned int[order];
    shape[0] = 1;
    shape[1] = 2;
    PeptideStateVector fpsv1(order, shape);
    fpsv1.tensor[{0, 0}] = 0.31;
    fpsv1.tensor[{0, 1}] = 0.71;
    PeptideStateVector bpsv1(order, shape);
    bpsv1.tensor[{0, 0}] = 0.32;
    bpsv1.tensor[{0, 1}] = 0.72;
    PeptideStateVector nbpsv1(order, shape);
    nbpsv1.tensor[{0, 0}] = 0.33;
    nbpsv1.tensor[{0, 1}] = 0.73;
    PeptideStateVector fpsv2(order, shape);
    fpsv2.tensor[{0, 0}] = 0.231;
    fpsv2.tensor[{0, 1}] = 0.271;
    PeptideStateVector bpsv2(order, shape);
    bpsv2.tensor[{0, 0}] = 0.232;
    bpsv2.tensor[{0, 1}] = 0.272;
    PeptideStateVector nbpsv2(order, shape);
    nbpsv2.tensor[{0, 0}] = 0.233;
    nbpsv2.tensor[{0, 1}] = 0.273;
    delete[] shape;
    unsigned int edmans = 0;
    double prob1 = 0.123456789;
    double prob2 = 0.987654321;
    ParameterFitter pf;
    bt.improve_fit(fpsv1, bpsv1, nbpsv1, edmans, prob1, &pf);
    bt.improve_fit(fpsv2, bpsv2, nbpsv2, edmans, prob2, &pf);
    BOOST_TEST(pf.get()
               == (0.71 * q * 0.33 / prob1 + 0.271 * q * 0.233 / prob2)
                          / (0.71 * 0.72 / prob1 + 0.271 * 0.272 / prob2));
}

BOOST_AUTO_TEST_SUITE_END()  // binomial_transition_suite
BOOST_AUTO_TEST_SUITE_END()  // step_suite
BOOST_AUTO_TEST_SUITE_END()  // hmm_suite

}  // namespace whatprot
