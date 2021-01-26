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

namespace fluoroseq {

namespace {
using boost::unit_test::tolerance;
const double TOL = 0.000000001;
}  // namespace

BOOST_AUTO_TEST_SUITE(hmm_suite)
BOOST_AUTO_TEST_SUITE(binomial_transition_suite)

BOOST_AUTO_TEST_CASE(constructor_test, *tolerance(TOL)) {
    double q = 0.2;
    int channel = -1;  // can be ignored for this test.
    BinomialTransition bt(q, channel);
    BOOST_TEST(bt.q == q);
}

BOOST_AUTO_TEST_CASE(reserve_zero_test, *tolerance(TOL)) {
    double q = 0.2;
    int channel = -1;  // can be ignored for this test.
    BinomialTransition bt(q, channel);
    bt.reserve(0);
    BOOST_TEST(bt.prob(0, 0) == 1.0);
}

BOOST_AUTO_TEST_CASE(reserve_one_test, *tolerance(TOL)) {
    double q = 0.2;
    double p = 0.8;
    int channel = -1;  // can be ignored for this test.
    BinomialTransition bt(q, channel);
    bt.reserve(1);
    BOOST_TEST(bt.prob(0, 0) == 1.0);
    BOOST_TEST(bt.prob(1, 0) == q);
    BOOST_TEST(bt.prob(1, 1) == p);
}

BOOST_AUTO_TEST_CASE(reserve_two_test, *tolerance(TOL)) {
    double q = 0.2;
    double p = 0.8;
    int channel = -1;  // can be ignored for this test.
    BinomialTransition bt(q, channel);
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
    BinomialTransition bt(q, channel);
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
    BinomialTransition bt(q, channel);
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

BOOST_AUTO_TEST_CASE(forward_trivial_test, *tolerance(TOL)) {
    double q = 0.05;
    double p = 0.95;
    int channel = 0;
    BinomialTransition bt(q, channel);
    bt.reserve(0);
    int order = 2;
    int* shape = new int[order];
    shape[0] = 1;
    shape[1] = 1;
    Tensor tsr(order, shape);
    delete[] shape;
    int* loc = new int[order];
    loc[0] = 0;
    loc[1] = 0;
    tsr[loc] = 1.0;  // loc is {0, 0}
    int edmans = 0;
    bt.forward(tsr, &edmans, &tsr);
    BOOST_TEST(tsr[loc] == 1.0);  // loc is {0, 0}
    delete[] loc;
}

BOOST_AUTO_TEST_CASE(forward_basic_transition_test, *tolerance(TOL)) {
    double q = 0.05;
    double p = 0.95;
    int channel = 0;
    BinomialTransition bt(q, channel);
    bt.reserve(1);
    int order = 2;
    int* shape = new int[order];
    shape[0] = 1;
    shape[1] = 2;
    Tensor tsr(order, shape);
    delete[] shape;
    int* loc = new int[order];
    loc[0] = 0;
    loc[1] = 0;
    tsr[loc] = 0.3;  // loc is {0, 0}
    loc[1] = 1;
    tsr[loc] = 0.7;  // loc is {0, 1}
    int edmans = 0;
    bt.forward(tsr, &edmans, &tsr);
    loc[0] = 0;
    loc[1] = 0;
    BOOST_TEST(tsr[loc] == 0.3 + 0.7 * q);  // loc is {0, 0}
    loc[1] = 1;
    BOOST_TEST(tsr[loc] == 0.7 * p);  // loc is {0, 1}
    delete[] loc;
}

BOOST_AUTO_TEST_CASE(forward_bigger_transition_test, *tolerance(TOL)) {
    double q = 0.05;
    double p = 0.95;
    int channel = 0;
    BinomialTransition bt(q, channel);
    bt.reserve(2);
    int order = 2;
    int* shape = new int[order];
    shape[0] = 1;
    shape[1] = 3;
    Tensor tsr(order, shape);
    delete[] shape;
    int* loc = new int[order];
    loc[0] = 0;
    loc[1] = 0;
    tsr[loc] = 0.2;  // loc is {0, 0}
    loc[1] = 1;
    tsr[loc] = 0.3;  // loc is {0, 1}
    loc[1] = 2;
    tsr[loc] = 0.7;  // loc is {0, 2}
    int edmans = 0;
    bt.forward(tsr, &edmans, &tsr);
    loc[0] = 0;
    loc[1] = 0;
    BOOST_TEST(tsr[loc] == 0.2 + 0.3 * q + 0.7 * q * q);  // loc is {0, 0}
    loc[1] = 1;
    BOOST_TEST(tsr[loc] == 0.3 * p + 0.7 * 2 * q * p);  // loc is {0, 1}
    loc[1] = 2;
    BOOST_TEST(tsr[loc] == 0.7 * p * p);  // loc is {0, 2}
    delete[] loc;
}

BOOST_AUTO_TEST_CASE(forward_multiple_edmans_test, *tolerance(TOL)) {
    double q = 0.05;
    double p = 0.95;
    int channel = 0;
    BinomialTransition bt(q, channel);
    bt.reserve(2);
    int order = 2;
    int* shape = new int[order];
    shape[0] = 3;
    shape[1] = 2;
    Tensor tsr(order, shape);
    delete[] shape;
    int* loc = new int[order];
    loc[0] = 0;
    loc[1] = 0;
    tsr[loc] = 0.2;  // loc is {0, 0}
    loc[1] = 1;
    tsr[loc] = 0.8;  // loc is {0, 1}
    loc[0] = 1;
    loc[1] = 0;
    tsr[loc] = 0.3;  // loc is {1, 0}
    loc[1] = 1;
    tsr[loc] = 0.7;  // loc is {1, 1}
    loc[0] = 2;
    loc[1] = 0;
    tsr[loc] = 0.4;  // loc is {2, 0}
    loc[1] = 1;
    tsr[loc] = 0.6;  // loc is {2, 1}
    int edmans = 2;
    bt.forward(tsr, &edmans, &tsr);
    loc[0] = 0;
    loc[1] = 0;
    BOOST_TEST(tsr[loc] == 0.2 + 0.8 * q);  // loc is {0, 0}
    loc[1] = 1;
    BOOST_TEST(tsr[loc] == 0.8 * p);  // loc is {0, 1}
    loc[0] = 1;
    loc[1] = 0;
    BOOST_TEST(tsr[loc] == 0.3 + 0.7 * q);  // loc is {1, 0}
    loc[1] = 1;
    BOOST_TEST(tsr[loc] == 0.7 * p);  // loc is {1, 1}
    loc[0] = 2;
    loc[1] = 0;
    BOOST_TEST(tsr[loc] == 0.4 + 0.6 * q);  // loc is {2, 0}
    loc[1] = 1;
    BOOST_TEST(tsr[loc] == 0.6 * p);  // loc is {2, 1}
    delete[] loc;
}

BOOST_AUTO_TEST_CASE(forward_other_dye_colors_test, *tolerance(TOL)) {
    double q = 0.05;
    double p = 0.95;
    int channel = 1;  // corresponds to 2nd dim of tensor
    BinomialTransition bt(q, channel);
    bt.reserve(1);
    int order = 4;
    int* shape = new int[order];
    shape[0] = 1;
    shape[1] = 2;
    shape[2] = 2;
    shape[3] = 2;
    Tensor tsr(order, shape);
    delete[] shape;
    int* loc = new int[order];
    loc[0] = 0;
    loc[1] = 0;
    loc[2] = 0;
    loc[3] = 0;
    tsr[loc] = 0.1;  // loc is {0, 0, 0, 0}
    loc[3] = 1;
    tsr[loc] = 0.2;  // loc is {0, 0, 0, 1}
    loc[2] = 1;
    loc[3] = 0;
    tsr[loc] = 0.3;  // loc is {0, 0, 1, 0}
    loc[3] = 1;
    tsr[loc] = 0.4;  // loc is {0, 0, 1, 1}
    loc[1] = 1;
    loc[2] = 0;
    loc[3] = 0;
    tsr[loc] = 0.5;  // loc is {0, 1, 0, 0}
    loc[3] = 1;
    tsr[loc] = 0.6;  // loc is {0, 1, 0, 1}
    loc[2] = 1;
    loc[3] = 0;
    tsr[loc] = 0.7;  // loc is {0, 1, 1, 0}
    loc[3] = 1;
    tsr[loc] = 0.8;  // loc is {0, 1, 1, 1}
    int edmans = 0;
    bt.forward(tsr, &edmans, &tsr);
    loc[0] = 0;
    loc[1] = 0;
    loc[2] = 0;
    loc[3] = 0;
    BOOST_TEST(tsr[loc] == 0.1 + 0.3 * q);  // loc is {0, 0, 0, 0}
    loc[3] = 1;
    BOOST_TEST(tsr[loc] == 0.2 + 0.4 * q);  // loc is {0, 0, 0, 1}
    loc[2] = 1;
    loc[3] = 0;
    BOOST_TEST(tsr[loc] == 0.3 * p);  // loc is {0, 0, 1, 0}
    loc[3] = 1;
    BOOST_TEST(tsr[loc] == 0.4 * p);  // loc is {0, 0, 1, 1}
    loc[1] = 1;
    loc[2] = 0;
    loc[3] = 0;
    BOOST_TEST(tsr[loc] == 0.5 + 0.7 * q);  // loc is {0, 1, 0, 0}
    loc[3] = 1;
    BOOST_TEST(tsr[loc] == 0.6 + 0.8 * q);  // loc is {0, 1, 0, 1}
    loc[2] = 1;
    loc[3] = 0;
    BOOST_TEST(tsr[loc] == 0.7 * p);  // loc is {0, 1, 1, 0}
    loc[3] = 1;
    BOOST_TEST(tsr[loc] == 0.8 * p);  // loc is {0, 1, 1, 1}
    delete[] loc;
}

BOOST_AUTO_TEST_CASE(backward_trivial_test, *tolerance(TOL)) {
    double q = 0.05;
    double p = 0.95;
    int channel = 0;
    BinomialTransition bt(q, channel);
    bt.reserve(0);
    int order = 2;
    int* shape = new int[order];
    shape[0] = 1;
    shape[1] = 1;
    Tensor tsr(order, shape);
    delete[] shape;
    int* loc = new int[order];
    loc[0] = 0;
    loc[1] = 0;
    tsr[loc] = 1.0;  // loc is {0, 0}
    int edmans = 0;
    bt.backward(tsr, &edmans, &tsr);
    BOOST_TEST(tsr[loc] == 1.0);  // loc is {0, 0}
    delete[] loc;
}

BOOST_AUTO_TEST_CASE(backward_basic_transition_test, *tolerance(TOL)) {
    double q = 0.05;
    double p = 0.95;
    int channel = 0;
    BinomialTransition bt(q, channel);
    bt.reserve(1);
    int order = 2;
    int* shape = new int[order];
    shape[0] = 1;
    shape[1] = 2;
    Tensor tsr(order, shape);
    delete[] shape;
    int* loc = new int[order];
    loc[0] = 0;
    loc[1] = 0;
    tsr[loc] = 0.3;  // loc is {0, 0}
    loc[1] = 1;
    tsr[loc] = 0.7;  // loc is {0, 1}
    int edmans = 0;
    bt.backward(tsr, &edmans, &tsr);
    loc[0] = 0;
    loc[1] = 0;
    BOOST_TEST(tsr[loc] == 0.3);  // loc is {0, 0}
    loc[1] = 1;
    BOOST_TEST(tsr[loc] == q * 0.3 + p * 0.7);  // loc is {0, 1}
    delete[] loc;
}

BOOST_AUTO_TEST_CASE(backward_bigger_transition_test, *tolerance(TOL)) {
    double q = 0.05;
    double p = 0.95;
    int channel = 0;
    BinomialTransition bt(q, channel);
    bt.reserve(2);
    int order = 2;
    int* shape = new int[order];
    shape[0] = 1;
    shape[1] = 3;
    Tensor tsr(order, shape);
    delete[] shape;
    int* loc = new int[order];
    loc[0] = 0;
    loc[1] = 0;
    tsr[loc] = 0.2;  // loc is {0, 0}
    loc[1] = 1;
    tsr[loc] = 0.3;  // loc is {0, 1}
    loc[1] = 2;
    tsr[loc] = 0.7;  // loc is {0, 2}
    int edmans = 0;
    bt.backward(tsr, &edmans, &tsr);
    loc[0] = 0;
    loc[1] = 0;
    BOOST_TEST(tsr[loc] == 0.2);  // loc is {0, 0}
    loc[1] = 1;
    BOOST_TEST(tsr[loc] == q * 0.2 + p * 0.3);  // loc is {0, 1}
    loc[1] = 2;
    // loc is {0, 2}
    BOOST_TEST(tsr[loc] == q * q * 0.2 + 2 * q * p * 0.3 + p * p * 0.7);
    delete[] loc;
}

BOOST_AUTO_TEST_CASE(backward_multiple_edmans_test, *tolerance(TOL)) {
    double q = 0.05;
    double p = 0.95;
    int channel = 0;
    BinomialTransition bt(q, channel);
    bt.reserve(2);
    int order = 2;
    int* shape = new int[order];
    shape[0] = 3;
    shape[1] = 2;
    Tensor tsr(order, shape);
    delete[] shape;
    int* loc = new int[order];
    loc[0] = 0;
    loc[1] = 0;
    tsr[loc] = 0.2;  // loc is {0, 0}
    loc[1] = 1;
    tsr[loc] = 0.8;  // loc is {0, 1}
    loc[0] = 1;
    loc[1] = 0;
    tsr[loc] = 0.3;  // loc is {1, 0}
    loc[1] = 1;
    tsr[loc] = 0.7;  // loc is {1, 1}
    loc[0] = 2;
    loc[1] = 0;
    tsr[loc] = 0.4;  // loc is {2, 0}
    loc[1] = 1;
    tsr[loc] = 0.6;  // loc is {2, 1}
    int edmans = 2;
    bt.backward(tsr, &edmans, &tsr);
    loc[0] = 0;
    loc[1] = 0;
    BOOST_TEST(tsr[loc] == 0.2);  // loc is {0, 0}
    loc[1] = 1;
    BOOST_TEST(tsr[loc] == q * 0.2 + p * 0.8);  // loc is {0, 1}
    loc[0] = 1;
    loc[1] = 0;
    BOOST_TEST(tsr[loc] == 0.3);  // loc is {1, 0}
    loc[1] = 1;
    BOOST_TEST(tsr[loc] == q * 0.3 + p * 0.7);  // loc is {1, 1}
    loc[0] = 2;
    loc[1] = 0;
    BOOST_TEST(tsr[loc] == 0.4);  // loc is {2, 0}
    loc[1] = 1;
    BOOST_TEST(tsr[loc] == q * 0.4 + p * 0.6);  // loc is {2, 1}
    delete[] loc;
}

BOOST_AUTO_TEST_CASE(backward_other_dye_colors_test, *tolerance(TOL)) {
    double q = 0.05;
    double p = 0.95;
    int channel = 1;  // corresponds to 2nd dim of tensor
    BinomialTransition bt(q, channel);
    bt.reserve(1);
    int order = 4;
    int* shape = new int[order];
    shape[0] = 1;
    shape[1] = 2;
    shape[2] = 2;
    shape[3] = 2;
    Tensor tsr(order, shape);
    delete[] shape;
    int* loc = new int[order];
    loc[0] = 0;
    loc[1] = 0;
    loc[2] = 0;
    loc[3] = 0;
    tsr[loc] = 0.1;  // loc is {0, 0, 0, 0}
    loc[3] = 1;
    tsr[loc] = 0.2;  // loc is {0, 0, 0, 1}
    loc[2] = 1;
    loc[3] = 0;
    tsr[loc] = 0.3;  // loc is {0, 0, 1, 0}
    loc[3] = 1;
    tsr[loc] = 0.4;  // loc is {0, 0, 1, 1}
    loc[1] = 1;
    loc[2] = 0;
    loc[3] = 0;
    tsr[loc] = 0.5;  // loc is {0, 1, 0, 0}
    loc[3] = 1;
    tsr[loc] = 0.6;  // loc is {0, 1, 0, 1}
    loc[2] = 1;
    loc[3] = 0;
    tsr[loc] = 0.7;  // loc is {0, 1, 1, 0}
    loc[3] = 1;
    tsr[loc] = 0.8;  // loc is {0, 1, 1, 1}
    int edmans = 0;
    bt.backward(tsr, &edmans, &tsr);
    loc[0] = 0;
    loc[1] = 0;
    loc[2] = 0;
    loc[3] = 0;
    BOOST_TEST(tsr[loc] == 0.1);  // loc is {0, 0, 0, 0}
    loc[3] = 1;
    BOOST_TEST(tsr[loc] == 0.2);  // loc is {0, 0, 0, 1}
    loc[2] = 1;
    loc[3] = 0;
    BOOST_TEST(tsr[loc] == q * 0.1 + p * 0.3);  // loc is {0, 0, 1, 0}
    loc[3] = 1;
    BOOST_TEST(tsr[loc] == q * 0.2 + p * 0.4);  // loc is {0, 0, 1, 1}
    loc[1] = 1;
    loc[2] = 0;
    loc[3] = 0;
    BOOST_TEST(tsr[loc] == 0.5);  // loc is {0, 1, 0, 0}
    loc[3] = 1;
    BOOST_TEST(tsr[loc] == 0.6);  // loc is {0, 1, 0, 1}
    loc[2] = 1;
    loc[3] = 0;
    BOOST_TEST(tsr[loc] == q * 0.5 + p * 0.7);  // loc is {0, 1, 1, 0}
    loc[3] = 1;
    BOOST_TEST(tsr[loc] == q * 0.6 + p * 0.8);  // loc is {0, 1, 1, 1}
    delete[] loc;
}

BOOST_AUTO_TEST_SUITE_END()  // binomial_transition_suite
BOOST_AUTO_TEST_SUITE_END()  // hmm_suite

}  // namespace fluoroseq
