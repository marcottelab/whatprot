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
#include "detach-transition.h"

namespace fluoroseq {

namespace {
using boost::unit_test::tolerance;
const double TOL = 0.000000001;
}  // namespace

BOOST_AUTO_TEST_SUITE(hmm_suite)
BOOST_AUTO_TEST_SUITE(detach_transition_suite)

BOOST_AUTO_TEST_CASE(constructor_test, *tolerance(TOL)) {
    double p_detach = 0.05;
    DetachTransition dt(p_detach);
    BOOST_CHECK_EQUAL(dt.p_detach, p_detach);
}

BOOST_AUTO_TEST_CASE(forward_trivial_test, *tolerance(TOL)) {
    double p_detach = 0.05;
    DetachTransition dt(p_detach);
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
    dt.forward(tsr, &edmans, &tsr);
    BOOST_TEST(tsr[loc] == 1.0);  // loc is {0, 0}
    delete[] loc;
}

BOOST_AUTO_TEST_CASE(forward_basic_test, *tolerance(TOL)) {
    double p_detach = 0.05;
    DetachTransition dt(p_detach);
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
    dt.forward(tsr, &edmans, &tsr);
    loc[0] = 0;
    loc[1] = 0;
    BOOST_TEST(tsr[loc] == 0.3 + 0.7 * p_detach);  // loc is {0, 0}
    loc[1] = 1;
    BOOST_TEST(tsr[loc] == 0.7 * (1 - p_detach));  // loc is {0, 1}
    delete[] loc;
}

BOOST_AUTO_TEST_CASE(forward_bigger_test, *tolerance(TOL)) {
    double p_detach = 0.05;
    DetachTransition dt(p_detach);
    int order = 2;
    int* shape = new int[order];
    shape[0] = 1;
    shape[1] = 3;
    Tensor tsr(order, shape);
    delete[] shape;
    int* loc = new int[order];
    loc[0] = 0;
    loc[1] = 0;
    tsr[loc] = 0.3;  // loc is {0, 0}
    loc[1] = 1;
    tsr[loc] = 0.6;  // loc is {0, 1}
    loc[1] = 2;
    tsr[loc] = 0.1;  // loc is {0, 2}
    int edmans = 0;
    dt.forward(tsr, &edmans, &tsr);
    loc[0] = 0;
    loc[1] = 0;
    BOOST_TEST(tsr[loc] == 0.3 + (0.6 + 0.1) * p_detach);  // loc is {0, 0}
    loc[1] = 1;
    BOOST_TEST(tsr[loc] == 0.6 * (1 - p_detach));  // loc is {0, 1}
    loc[1] = 2;
    BOOST_TEST(tsr[loc] == 0.1 * (1 - p_detach));  // loc is {0, 2}
    delete[] loc;
}

BOOST_AUTO_TEST_CASE(forward_multiple_edmans_test, *tolerance(TOL)) {
    double p_detach = 0.05;
    DetachTransition dt(p_detach);
    int order = 2;
    int* shape = new int[order];
    shape[0] = 3;
    shape[1] = 2;
    Tensor tsr(order, shape);
    delete[] shape;
    int* loc = new int[order];
    loc[0] = 0;
    loc[1] = 0;
    tsr[loc] = 0.1;  // loc is {0, 0}
    loc[1] = 1;
    tsr[loc] = 0.2;  // loc is {0, 1}
    loc[0] = 1;
    loc[1] = 0;
    tsr[loc] = 0.3;  // loc is {1, 0}
    loc[1] = 1;
    tsr[loc] = 0.4;  // loc is {1, 1}
    loc[0] = 2;
    loc[1] = 0;
    tsr[loc] = 0.5;  // loc is {2, 0}
    loc[1] = 1;
    tsr[loc] = 0.6;  // loc is {2, 1}
    int edmans = 2;
    dt.forward(tsr, &edmans, &tsr);
    // Just testing the ones with at least one lit amino acid here. See below
    // for other tests.
    loc[0] = 0;
    loc[1] = 1;
    BOOST_TEST(tsr[loc] == 0.2 * (1 - p_detach));  // loc is {0, 1}
    loc[0] = 1;
    BOOST_TEST(tsr[loc] == 0.4 * (1 - p_detach));  // loc is {1, 1}
    loc[0] = 2;
    BOOST_TEST(tsr[loc] == 0.6 * (1 - p_detach));  // loc is {2, 1}
    // Distribution between empty states is of no importance. We only care about
    // the sum.
    double sum_empties = 0.0;
    loc[0] = 0;
    loc[1] = 0;
    sum_empties += tsr[loc];  // loc is {0, 0}
    loc[0] = 1;
    sum_empties += tsr[loc];  // loc is {1, 0};
    loc[0] = 2;
    sum_empties += tsr[loc];  // loc is {2, 0};
    BOOST_TEST(sum_empties == 0.1 + 0.3 + 0.5 + (0.2 + 0.4 + 0.6) * p_detach);
    delete[] loc;
}

BOOST_AUTO_TEST_CASE(forward_multiple_dye_colors_test, *tolerance(TOL)) {
    double p_detach = 0.05;
    DetachTransition dt(p_detach);
    int order = 3;
    int* shape = new int[order];
    shape[0] = 1;
    shape[1] = 2;
    shape[2] = 2;
    Tensor tsr(order, shape);
    delete[] shape;
    int* loc = new int[order];
    loc[0] = 0;
    loc[1] = 0;
    loc[2] = 0;
    tsr[loc] = 0.1;  // loc is {0, 0, 0}
    loc[2] = 1;
    tsr[loc] = 0.2;  // loc is {0, 0, 1}
    loc[1] = 1;
    loc[2] = 0;
    tsr[loc] = 0.3;  // loc is {0, 1, 0}
    loc[2] = 1;
    tsr[loc] = 0.4;  // loc is {0, 1, 1}
    int edmans = 0;
    dt.forward(tsr, &edmans, &tsr);
    loc[0] = 0;
    loc[1] = 0;
    loc[2] = 0;
    // loc is {0, 0, 0}
    BOOST_TEST(tsr[loc] == 0.1 + (0.2 + 0.3 + 0.4) * p_detach);
    loc[2] = 1;
    BOOST_TEST(tsr[loc] == 0.2 * (1 - p_detach));  // loc is {0, 0, 1}
    loc[1] = 1;
    loc[2] = 0;
    BOOST_TEST(tsr[loc] == 0.3 * (1 - p_detach));  // loc is {0, 1, 0}
    loc[2] = 1;
    BOOST_TEST(tsr[loc] == 0.4 * (1 - p_detach));  // loc is {0, 1, 1}
    delete[] loc;
}

BOOST_AUTO_TEST_CASE(backward_trivial_test, *tolerance(TOL)) {
    double p_detach = 0.05;
    DetachTransition dt(p_detach);
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
    dt.backward(tsr, &edmans, &tsr);
    BOOST_TEST(tsr[loc] == 1.0);  // loc is {0, 0}
    delete[] loc;
}

BOOST_AUTO_TEST_CASE(backward_basic_test, *tolerance(TOL)) {
    double p_detach = 0.05;
    DetachTransition dt(p_detach);
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
    dt.backward(tsr, &edmans, &tsr);
    loc[0] = 0;
    loc[1] = 0;
    BOOST_TEST(tsr[loc] == 0.3);  // loc is {0, 0}
    loc[1] = 1;
    // loc is {0, 1}
    BOOST_TEST(tsr[loc] == p_detach * 0.3 + (1 - p_detach) * 0.7);
    delete[] loc;
}

BOOST_AUTO_TEST_CASE(backward_bigger_test, *tolerance(TOL)) {
    double p_detach = 0.05;
    DetachTransition dt(p_detach);
    int order = 2;
    int* shape = new int[order];
    shape[0] = 1;
    shape[1] = 3;
    Tensor tsr(order, shape);
    delete[] shape;
    int* loc = new int[order];
    loc[0] = 0;
    loc[1] = 0;
    tsr[loc] = 0.3;  // loc is {0, 0}
    loc[1] = 1;
    tsr[loc] = 0.6;  // loc is {0, 1}
    loc[1] = 2;
    tsr[loc] = 0.1;  // loc is {0, 2}
    int edmans = 0;
    dt.backward(tsr, &edmans, &tsr);
    loc[0] = 0;
    loc[1] = 0;
    BOOST_TEST(tsr[loc] == 0.3);  // loc is {0, 0}
    loc[1] = 1;
    // loc is {0, 1}
    BOOST_TEST(tsr[loc] == p_detach * 0.3 + (1 - p_detach) * 0.6);
    loc[1] = 2;
    // loc is {0, 2}
    BOOST_TEST(tsr[loc] == p_detach * 0.3 + (1 - p_detach) * 0.1);
    delete[] loc;
}

BOOST_AUTO_TEST_CASE(backward_multiple_edmans_test, *tolerance(TOL)) {
    double p_detach = 0.05;
    DetachTransition dt(p_detach);
    int order = 2;
    int* shape = new int[order];
    shape[0] = 3;
    shape[1] = 2;
    Tensor tsr(order, shape);
    delete[] shape;
    int* loc = new int[order];
    loc[0] = 0;
    loc[1] = 0;
    tsr[loc] = 0.88;  // loc is {0, 0}
    loc[1] = 1;
    tsr[loc] = 0.2;  // loc is {0, 1}
    loc[0] = 1;
    loc[1] = 0;
    tsr[loc] = 0.88;  // loc is {1, 0}
    loc[1] = 1;
    tsr[loc] = 0.4;  // loc is {1, 1}
    loc[0] = 2;
    loc[1] = 0;
    tsr[loc] = 0.88;  // loc is {2, 0}
    loc[1] = 1;
    tsr[loc] = 0.6;  // loc is {2, 1}
    int edmans = 2;
    dt.backward(tsr, &edmans, &tsr);
    // Just testing the ones with at least one lit amino acid here. See below
    // for other tests.
    loc[0] = 0;
    loc[1] = 0;
    BOOST_TEST(tsr[loc] == 0.88);  // loc is {0, 0}
    loc[1] = 1;
    // loc is {0, 1}
    BOOST_TEST(tsr[loc] == p_detach * 0.88 + (1 - p_detach) * 0.2);
    loc[0] = 1;
    loc[1] = 0;
    BOOST_TEST(tsr[loc] == 0.88);  // loc is {1, 0}
    loc[1] = 1;
    // loc is {1, 1}
    BOOST_TEST(tsr[loc] == p_detach * 0.88 + (1 - p_detach) * 0.4);
    loc[0] = 2;
    loc[1] = 0;
    BOOST_TEST(tsr[loc] == 0.88);  // loc is {2, 0}
    loc[1] = 1;
    // loc is {2, 1}
    BOOST_TEST(tsr[loc] == p_detach * 0.88 + (1 - p_detach) * 0.6);
    delete[] loc;
}

BOOST_AUTO_TEST_CASE(backward_multiple_dye_colors_test, *tolerance(TOL)) {
    double p_detach = 0.05;
    DetachTransition dt(p_detach);
    int order = 3;
    int* shape = new int[order];
    shape[0] = 1;
    shape[1] = 2;
    shape[2] = 2;
    Tensor tsr(order, shape);
    delete[] shape;
    int* loc = new int[order];
    loc[0] = 0;
    loc[1] = 0;
    loc[2] = 0;
    tsr[loc] = 0.1;  // loc is {0, 0, 0}
    loc[2] = 1;
    tsr[loc] = 0.2;  // loc is {0, 0, 1}
    loc[1] = 1;
    loc[2] = 0;
    tsr[loc] = 0.3;  // loc is {0, 1, 0}
    loc[2] = 1;
    tsr[loc] = 0.4;  // loc is {0, 1, 1}
    int edmans = 0;
    dt.backward(tsr, &edmans, &tsr);
    loc[0] = 0;
    loc[1] = 0;
    loc[2] = 0;
    BOOST_TEST(tsr[loc] == 0.1);  // loc is {0, 0, 0}
    loc[2] = 1;
    // loc is {0, 0, 1}
    BOOST_TEST(tsr[loc] == p_detach * 0.1 + (1 - p_detach) * 0.2);
    loc[1] = 1;
    loc[2] = 0;
    // loc is {0, 1, 0}
    BOOST_TEST(tsr[loc] == p_detach * 0.1 + (1 - p_detach) * 0.3);
    loc[2] = 1;
    // loc is {0, 1, 1}
    BOOST_TEST(tsr[loc] == p_detach * 0.1 + (1 - p_detach) * 0.4);
    delete[] loc;
}

BOOST_AUTO_TEST_SUITE_END()  // detach_transition_suite
BOOST_AUTO_TEST_SUITE_END()  // hmm_suite

}  // namespace fluoroseq
