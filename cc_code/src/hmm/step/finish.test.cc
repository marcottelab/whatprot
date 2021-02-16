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
#include "finish.h"

// Local project headers:
#include "tensor/tensor.h"

namespace whatprot {

namespace {
using boost::unit_test::tolerance;
const double TOL = 0.000000001;
}  // namespace

BOOST_AUTO_TEST_SUITE(hmm_suite)
BOOST_AUTO_TEST_SUITE(step_suite)
BOOST_AUTO_TEST_SUITE(finish_suite)

BOOST_AUTO_TEST_CASE(forward_test, *tolerance(TOL)) {
    Finish finish;
    int order = 3;
    int* shape = new int[order];
    shape[0] = 2;
    shape[1] = 2;
    shape[2] = 2;
    Tensor tsr1(order, shape);
    Tensor tsr2(order, shape);
    delete[] shape;
    int* loc = new int[order];
    loc[0] = 0;
    loc[1] = 0;
    loc[2] = 0;
    tsr1[loc] = 0.00;  // loc is {0, 0, 0}
    loc[2] = 1;
    tsr1[loc] = 0.01;  // loc is {0, 0, 1}
    loc[1] = 1;
    loc[2] = 0;
    tsr1[loc] = 0.10;  // loc is {0, 1, 0}
    loc[2] = 1;
    tsr1[loc] = 0.11;  // loc is {0, 1, 1}
    int edmans;  // should be ignored.
    finish.forward(tsr1, &edmans, &tsr2);
    loc[1] = 0;
    loc[2] = 0;
    BOOST_TEST(tsr2[loc] == 0.00);  // loc is {0, 0, 0}
    loc[2] = 1;
    BOOST_TEST(tsr2[loc] == 0.01);  // loc is {0, 0, 1}
    loc[1] = 1;
    loc[2] = 0;
    BOOST_TEST(tsr2[loc] == 0.10);  // loc is {0, 1, 0}
    loc[2] = 1;
    BOOST_TEST(tsr2[loc] == 0.11);  // loc is {0, 1, 1}
    delete[] loc;
}

BOOST_AUTO_TEST_CASE(backward_test, *tolerance(TOL)) {
    Finish finish;
    int order = 3;
    int* shape = new int[order];
    shape[0] = 2;
    shape[1] = 2;
    shape[2] = 2;
    Tensor tsr1(order, shape);
    Tensor tsr2(order, shape);
    delete[] shape;
    int* loc = new int[order];
    loc[0] = 0;
    loc[1] = 0;
    loc[2] = 0;
    tsr1[loc] = 0.0;  // loc is {0, 0, 0}
    loc[2] = 1;
    tsr1[loc] = 0.0;  // loc is {0, 0, 1}
    loc[1] = 1;
    loc[2] = 0;
    tsr1[loc] = 0.0;  // loc is {0, 1, 0}
    loc[2] = 1;
    tsr1[loc] = 0.0;  // loc is {0, 1, 1}
    int edmans;  // should be ignored.
    finish.backward(tsr1, &edmans, &tsr2);
    loc[1] = 0;
    loc[2] = 0;
    BOOST_TEST(tsr2[loc] == 1.0);  // loc is {0, 0, 0}
    loc[2] = 1;
    BOOST_TEST(tsr2[loc] == 1.0);  // loc is {0, 0, 1}
    loc[1] = 1;
    loc[2] = 0;
    BOOST_TEST(tsr2[loc] == 1.0);  // loc is {0, 1, 0}
    loc[2] = 1;
    BOOST_TEST(tsr2[loc] == 1.0);  // loc is {0, 1, 1}
    delete[] loc;
}

BOOST_AUTO_TEST_SUITE_END()  // finish_suite
BOOST_AUTO_TEST_SUITE_END()  // step_suite
BOOST_AUTO_TEST_SUITE_END()  // hmm_suite

}  // namespace whatprot
